"""Deterministic, evidence-linked interpretation. No Promethease dependency."""
from collections import Counter
from datetime import datetime, timezone
from importlib.resources import files
import json
from pathlib import Path
import re
from urllib.parse import unquote
from . import __version__
from .parser import parse_genome
from .references import connect_read

CATEGORIES = ['Review with a professional','Conflicting evidence','Uncertain significance','Other association','Benign / likely benign','Unclassified']

def display(value):
    return unquote(value).replace('_',' ')

def stars(review):
    if 'practice_guideline' in review:return 4
    if 'reviewed_by_expert_panel' in review:return 3
    if 'multiple_submitters' in review and 'no_conflicts' in review:return 2
    if 'criteria_provided' in review:return 1
    return 0

def classification_category(sig):
    normalized=sig.lower().replace(' ','_').replace(',_low_penetrance','')
    terms={term.strip('_') for term in re.split(r'[/|,;]',normalized)}
    if any('conflict' in s for s in terms):return CATEGORIES[1]
    if terms & {'pathogenic','likely_pathogenic','pathogenic_low_penetrance','likely_pathogenic_low_penetrance'}:
        if terms & {'benign','likely_benign','uncertain_significance'}:return CATEGORIES[1]
        return CATEGORIES[0]
    if any('uncertain' in t for t in terms):return CATEGORIES[2]
    if terms and terms <= {'benign','likely_benign'}:return CATEGORIES[4]
    if any(any(k in t for k in ['risk','association','drug_response','protective','affects']) for t in terms):return CATEGORIES[3]
    return CATEGORIES[5]

def explain_clinvar(r,gt):
    cat=classification_category(r['classification'])
    if r['conflicts'] and cat==CATEGORIES[0]:cat=CATEGORIES[1]
    if cat==CATEGORIES[0]:
        text='The file contains an allele with a pathogenic or likely pathogenic ClinVar classification. A raw-array call is not a diagnosis. Independent clinical confirmation and review of the exact condition, inheritance, penetrance, and assay reliability are needed before a health decision.'
    elif cat==CATEGORIES[1]:
        text='ClinVar records conflicting interpretations. This app does not resolve that disagreement or classify the result as disease-causing.'
    elif cat==CATEGORIES[2]:
        text='The clinical meaning of this allele is uncertain. An uncertain classification does not establish disease or justify a treatment change.'
    elif cat==CATEGORIES[4]:
        text='ClinVar classifies this allele as benign or likely benign in its submitted clinical context. This is not an all-clear for the gene or for conditions not covered by this record.'
    elif cat==CATEGORIES[3]:
        text='This is an association, susceptibility, protective, or drug-response annotation. It is not a diagnosis, a complete pharmacogenetic result, or an absolute-risk estimate.'
    else:
        text='The annotation is retained for inspection. No supported clinical interpretation is assigned by this app.'
    return cat,text

def trait_results(g):
    data=json.loads(files('genetics_app').joinpath('data/traits.json').read_text())
    rows=[]
    for r in data['rules']:
        pos=r['positions'][g.build]
        status,gt,ids=g.locus(r['chrom'],pos)
        # Require identity and build/coordinate agreement; no silent alias conversion.
        if status=='called' and r['rsid'] not in ids:status='marker_identity_unverified'
        if status=='called' and not any(row.gt!='--' for row in g.by_id.get(r['rsid'],[])):status='no_call'
        if status=='called' and (len(gt)!=r['ploidy'] or not set(gt)<=set(r['alleles'])):status='unsupported_alleles_or_ploidy'
        explanation=r['explanations'].get(gt,'') if status=='called' else 'No interpretation: '+status.replace('_',' ')+'.'
        rows.append(dict(id=r['id'],title=r['title'],gene=r['gene'],rsid=r['rsid'],chrom=r['chrom'],pos=pos,genotype=gt,status=status,explanation=explanation,limits=r['limits'],evidence=r['evidence'],sources=r['sources']))
    return rows,data['version']

def clinvar_results(g,path):
    db,meta=connect_read(path)
    try:
        if meta.get('kind')!='ClinVar' or meta.get('build')!=g.build:raise ValueError('ClinVar reference build does not match the input.')
        db.execute('create temp table target(chrom text,pos integer,primary key(chrom,pos))')
        db.executemany('insert into target values(?,?)',g.by_position)
        counts=Counter();findings=[];exclusions=[]
        for row in db.execute('select v.* from target t join variants v indexed by locus on t.chrom=v.chrom and t.pos=v.pos'):
            r=dict(row);counts['position_overlap_records']+=1
            status,gt,ids=g.locus(r['chrom'],r['pos'])
            if status=='called':
                if len(r['ref'])!=1 or len(r['alt'])!=1 or not set(r['ref']+r['alt'])<=set('ACGT'):status='non_biallelic_snv_unresolved'
                elif not set(gt)<=set(r['ref']+r['alt']):status='allele_mismatch'
                elif r['alt'] not in gt:status='reference_only'
                # ORIGIN is bit-coded. Somatic-only annotations cannot become germline findings.
                elif r['origin'].isdigit() and int(r['origin']) & 2 and not int(r['origin']) & 1:status='somatic_only_unresolved'
                elif not r['classification']:status='no_germline_classification'
                else:status='alt_present'
            counts[status]+=1
            if status=='alt_present':
                cat,text=explain_clinvar(r,gt)
                findings.append(dict(marker_ids=ids,rsid=('rs'+r['rs']) if r['rs'].isdigit() else '',chrom=r['chrom'],pos=r['pos'],genotype=gt,ref=r['ref'],alt=r['alt'],alt_copies=gt.count(r['alt']),gene=display(r['gene']),classification=display(r['classification']),review=display(r['review']),stars=stars(r['review']),conditions=display(r['conditions']),category=cat,explanation=text,conflicts=display(r['conflicts']),variation_id=r['id'],source='https://www.ncbi.nlm.nih.gov/clinvar/variation/'+r['id']+'/' if r['id'].isdigit() else 'https://www.ncbi.nlm.nih.gov/clinvar/',context='Variant-level aggregate; classification is not mapped separately to each listed condition.'))
            elif status!='reference_only':
                exclusions.append(dict(marker_ids=ids,chrom=r['chrom'],pos=r['pos'],genotype=gt,variation_id=r['id'],reason=status))
        findings.sort(key=lambda r:(CATEGORIES.index(r['category']),-r['stars'],r['chrom'],r['pos']))
        counts['findings']=len(findings)
        return findings,exclusions,dict(counts),meta
    finally:db.close()


def sequence_clinvar_results(knowledge):
    """Clinical context for independently sequence-resolved indels only."""
    findings=[];seen=set()
    for finding in knowledge['findings']:
        resolved=finding.get('sequence_resolution')
        if not resolved or resolved.get('coding_status')!='polarity_invariant_heterozygote':continue
        r=resolved['definition']['clinical_record'];gt=resolved['sequence_call']
        count=gt.count(r['inserted'])
        if not count or not r['classification']:continue
        if r['origin'].isdigit() and int(r['origin'])&2 and not int(r['origin'])&1:continue
        signature=(r['chrom'],r['pos'],r['ref'],r['alt'],r['id'])
        if signature in seen:continue
        seen.add(signature);cat,text=explain_clinvar(r,gt)
        findings.append(dict(marker_ids=finding['raw_ids'],rsid=finding['id'],chrom=r['chrom'],pos=r['pos'],
            genotype=' / '.join(a or 'deletion' for a in gt),ref=r['ref'],alt=r['alt'],alt_copies=count,
            gene=display(r['gene']),classification=display(r['classification']),review=display(r['review']),stars=stars(r['review']),
            conditions=display(r['conditions']),category=cat,explanation=text,conflicts=display(r['conflicts']),variation_id=r['id'],
            source='https://www.ncbi.nlm.nih.gov/clinvar/variation/'+r['id']+'/',
            context='Allele-specific sequence resolution from an explicitly identified I/D assay; variant-level aggregate, not a diagnosis.',
            sequence_resolution=resolved['definition']['method'],audit_flags=list(finding['flags'])))
    return findings

def gwas_results(g,path,limit=2000):
    db,meta=connect_read(path)
    try:
        if meta.get('kind')!='GWAS Catalog':raise ValueError('Wrong GWAS reference database.')
        db.execute('create temp table target(rs text primary key)')
        ids=[(i,) for i,rr in g.by_id.items() if re.fullmatch(r'rs\d+',i) and i not in g.bad_ids and g.locus(rr[0].chrom,rr[0].pos)[0]=='called']
        db.executemany('insert into target values(?)',ids)
        query='from target t join associations a indexed by marker on t.rs=a.rs'
        total=db.execute('select count(*) '+query).fetchone()[0]
        rows=[dict(r) for r in db.execute('select a.* '+query+' order by a.p,a.rs,a.pmid,a.trait limit ?', (limit,))]
        for r in rows:
            r['source']='https://pubmed.ncbi.nlm.nih.gov/'+r['pmid']+'/' if r['pmid'].isdigit() else 'https://www.ebi.ac.uk/gwas/'
            r['interpretation']='Literature lead only. The effect allele, strand, phenotype units, ancestry applicability and LD have not been harmonized; no personal risk direction is inferred.'
        return rows,dict(total_overlap_records=total,displayed_records=len(rows),limit=limit,truncated=total>limit),meta
    finally:db.close()

def analyze(data,build='auto',strand='forward',clinvar=None,gwas=None,traits_only=False):
    g=parse_genome(data,build,strand)
    if not clinvar and not traits_only:raise ValueError('Prepare a matching ClinVar reference, or explicitly select traits-only mode.')
    traits,version=trait_results(g)
    report=dict(app_version=__version__,generated_utc=datetime.now(timezone.utc).isoformat(),build=g.build,build_source=g.build_source,input_sha256=g.sha256,qc=g.qc,traits=traits,traits_version=version,clinvar=[],exclusions=[],clinvar_counts={},gwas=[],gwas_counts={},references=[],mode='traits only' if not clinvar else 'ClinVar and curated traits',limits=[
        'Research and education, not a diagnosis. Raw-array rare calls can be wrong; important findings require independent clinical confirmation.',
        'Only forward-strand GRCh37/GRCh38 exports and exact allele-compatible biallelic SNVs are interpreted. Missing, conflicting, indel and unsupported calls stay unresolved.',
        'No finding does not rule out disease. Unmeasured variants, CNVs, structural variants, repeats and comprehensive pharmacogenetic diplotypes are not evaluated.',
        'ClinVar stars describe review of a classification, not accuracy of your raw genotype or personal disease probability.',
        'This app has no SNPedia/Promethease dependency, genoset engine, ancestry inference, imputation or calibrated polygenic risk in the report path.'
    ])
    if clinvar:
        report['clinvar'],report['exclusions'],report['clinvar_counts'],meta=clinvar_results(g,clinvar);report['references'].append(meta)
    if gwas:
        report['gwas'],report['gwas_counts'],meta=gwas_results(g,gwas);report['references'].append(meta)
    return report
