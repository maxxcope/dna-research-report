"""Complete local orchestration: generate, audit, compute, and explain."""
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile
import time
from .parser import parse_genome
from .engine import analyze, sequence_clinvar_results
from .knowledge import interpret
from . import extras, scoring
from .references import digest


def load_config(path):
    path=Path(path).resolve();config=json.loads(path.read_text())
    if config.get('schema')!=1:raise ValueError('Unsupported pipeline profile schema.')
    # Path values are local configuration, never supplied by an uploaded genome.
    keys={'clusters','badalleles','prefix','labels','plink','tree','chain','directory','java','beagle','score','file','map','path','knowledge','clinvar','gwas'}
    def expand(value,key=''):
        if isinstance(value,dict):return {k:expand(v,k) for k,v in value.items()}
        if isinstance(value,list):return [expand(v) for v in value]
        if key in keys and isinstance(value,str) and not value.startswith(('https://','http://')):
            # HaploGrep's tree ID isn't a filename.
            if key=='tree' and '@' in value:return value
            p=Path(value).expanduser();return str(p if p.is_absolute() else path.parent/p)
        return value
    return expand(config)


def full_report(data,profile,build='auto',progress=lambda _:None,work=None,events=lambda _:None):
    config=load_config(profile);genome=parse_genome(data,build);stages=[];started=time.monotonic()
    extras.reset_commands()
    total_stages = 10 + max(1, len(config.get('ancestry', [])))
    def emit(kind, name='', **values):
        events(dict(kind=kind, name=name, completed=len(stages), total=total_stages, **values))
    def detail(message):
        progress(message)
        events(dict(kind='detail', message=message))
    progress('Verifying public reference files')
    for index, reference in enumerate(config.get('reference_files',[]), 1):
        if digest(reference['file'])!=reference['sha256']:raise ValueError('Reference checksum changed: '+Path(reference['file']).name+'. Prepare a new profile before analysis.')
        events(dict(kind='references', current=index, total=len(config.get('reference_files', []))))
    def stage(name,fn):
        progress(name);emit('stage_started', name);t=time.monotonic()
        try:
            value=fn();status=value.get('status','completed') if isinstance(value,dict) else 'completed'
            stages.append(dict(name=name,status=status,seconds=round(time.monotonic()-t,2)))
            emit('stage_finished', name, status=status)
            return value
        except Exception as e:
            stages.append(dict(name=name,status='failed',seconds=round(time.monotonic()-t,2)))
            emit('stage_finished', name, status='failed')
            return dict(status='failed',explanation=str(e))
    def absent(name):
        emit('stage_started', name)
        stages.append(dict(name=name,status='reference_not_configured'))
        emit('stage_finished', name, status='reference_not_configured')
        return dict(status='reference_not_configured',explanation=f'{name} needs its public reference in the local pipeline profile.')
    flagged = set()
    def quality_check():
        nonlocal flagged
        result, flagged = extras.quality(genome, config['quality'])
        return result
    quality = stage('Array quality', quality_check) if config.get('quality') else absent('Array quality')
    cv=config.get('clinvar');gw=config.get('gwas')
    # Full report must have both independent annotation channels configured.
    if not cv or not config.get('knowledge'):raise ValueError('Full pipeline requires ClinVar and SNPedia references; run reference setup first.')
    report=stage('Clinical annotations and measured traits',lambda:analyze(data,build,clinvar=cv,gwas=gw))
    if report.get('status')=='failed':raise ValueError(report['explanation'])
    knowledge=stage('SNPedia findings and multi-marker rules',lambda:interpret(genome,config['knowledge'],flagged))
    if knowledge.get('status')=='failed':raise ValueError(knowledge['explanation'])
    report['knowledge']=knowledge;report['references'].append(knowledge['reference']);report['research']={};report['mode']='Full local research pipeline'
    sequence_annotations=sequence_clinvar_results(knowledge)
    report['clinvar'].extend(sequence_annotations)
    report['clinvar_counts']['sequence_resolved_indel_annotations']=len(sequence_annotations)
    report['clinvar_counts']['findings']=len(report['clinvar'])
    report['research']['array_quality']=quality
    for trait in report['traits']:
        trait['audit_flags']=['published_array_quality_flag'] if (trait['chrom'],trait['pos']) in flagged else []
    # Audit findings after generation against independent contemporary annotations.
    bypos={}
    for r in report['clinvar']:
        bypos.setdefault((r['chrom'],r['pos']),[]).append(r)
        r['audit_flags']=sorted(set(r.get('audit_flags',[])+(['published_array_quality_flag'] if (r['chrom'],r['pos']) in flagged else [])))
    audit=[]
    for finding in knowledge['findings']:
        annotations=[r for r in bypos.get((finding.get('chrom'),finding.get('pos')),[]) if r['rsid']==finding['id']]
        if finding.get('sequence_resolution'):
            vid=finding['sequence_resolution']['definition']['clinical_record']['id']
            annotations=[r for r in sequence_annotations if r['variation_id']==vid]
        finding['clinical_crosscheck']=[dict(id=r['variation_id'],classification=r['classification'],review=r['review'],source=r['source']) for r in annotations]
        if any(r['category']=='Review with a professional' for r in annotations) and finding['repute'].lower()=='good':
            finding['flags'].append('source_interpretations_need_context_review')
        if finding['flags']:audit.append(dict(id=finding['id'],flags=finding['flags']))
    report['finding_audit']=dict(flagged_findings=audit,unresolved_variants=knowledge['unresolved'],genosets=knowledge['genoset_audit'],explanation='Findings are checked for coordinate/strand compatibility, missing or conflicting calls, published array-quality flags and overlap with ClinVar. Differences between sources are retained for review, not silently resolved.')
    def calculate(workdir):
        research=report['research']
        research['ancestry']=[]
        for panel in config.get('ancestry',[]):
            research['ancestry'].append(stage('Ancestry: '+panel['name'],lambda panel=panel:extras.ancestry(genome,panel,workdir/('ancestry_'+str(len(research['ancestry']))),flagged)))
        if not research['ancestry']:research['ancestry']=[absent('Ancestry')]
        research['frequencies']=stage('Reference genotype frequencies',lambda:extras.frequencies(genome,config['frequencies'],knowledge['findings']) if genome.build=='37' else dict(status='unsupported_build',explanation='This reference panel uses GRCh37.')) if config.get('frequencies') else absent('Reference genotype frequencies')
        research['paternal_lineage']=stage('Paternal tree comparison',lambda:extras.y_lineage(genome,config['paternal'])) if config.get('paternal') else absent('Paternal lineage')
        research['maternal_lineage']=stage('Maternal lineage: HaploGrep',lambda:extras.mito(genome,config['maternal'],workdir/'mito')) if config.get('maternal') else absent('Maternal lineage')
        research['score_coverage']=stage('Published score coverage',lambda:scoring.direct_scores(genome,config.get('scores',[]),flagged))
        research['imputed_score']=stage('Local imputation and score',lambda:scoring.imputed_score(genome,config['imputation'],workdir/'imputation',flagged,detail)) if config.get('imputation') else absent('Imputed score')
        research['phase']=stage('Targeted statistical phasing',lambda:scoring.phase(genome,config['phase'],workdir/'phase',flagged)) if config.get('phase') else absent('Statistical phase')
    if work:
        workdir=Path(work).resolve();workdir.mkdir(parents=True,exist_ok=True)
        if any(workdir.iterdir()):raise ValueError('Use an empty work directory for each input; cached personal outputs cannot be reused.')
        calculate(workdir)
    else:
        with tempfile.TemporaryDirectory(prefix='genetics-run-') as tmp:calculate(Path(tmp))
    from .conclusions import make_report
    emit('stage_started', 'Finding audit and report assembly')
    report['conclusions']=make_report(report,genome)
    stages.append(dict(name='Finding audit and report assembly',status='completed'))
    report['pipeline']=dict(stages=stages,seconds=round(time.monotonic()-started,2),profile_sha256=digest(profile),status='completed' if all(s['status']=='completed' for s in stages) else 'completed_with_unresolved_stages')
    report['pipeline']['reference_files']=[dict(name=Path(r['file']).name,sha256=r['sha256']) for r in config.get('reference_files',[])]
    report['limits']=[x for x in report['limits'] if not x.startswith('This app has no SNPedia')]
    report['limits']=[x for x in report['limits'] if not x.startswith('Only forward-strand')]
    report['limits'].append('Matching supports compatible SNVs and simple heterozygous I/D sequence mappings that are unchanged by reversing provider codes. Homozygous I/D conversions without provider-specific allele definitions remain unresolved; literal source I/D pages do not establish sequence identity.')
    report['limits'] += ['SNPedia summaries, source labels and literature are reproduced with attribution under CC BY-NC-SA 3.0 US. Magnitude is subjective source interest, not disease probability.','Each research section reports its own coverage, reference panel and status. Statistical phase and imputed calls are labelled separately from measured genotypes.','Source agreement is an audit of interpretations. Missing rules, missing calls and unsupported variant representations remain visible in the audit.']
    emit('stage_finished', 'Finding audit and report assembly', status='completed')
    emit('complete', 'Analysis finished', status=report['pipeline']['status'])
    progress('Full report assembled')
    return report
