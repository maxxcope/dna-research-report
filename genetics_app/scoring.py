"""Published-score coverage, local Beagle inference and explicit phase evidence."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import hashlib
import json
from pathlib import Path
from .extras import run_command, cancel_commands


def score_rows(path):
    with gzip.open(path,'rt') as f:
        return list(csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t'))


def direct_scores(genome,definitions,flagged):
    results=[]
    for definition in definitions:
        if definition.get('build','37')!=genome.build:
            results.append(dict(id=definition['id'],status='unsupported_build'));continue
        rows=score_rows(definition['path']);components=[];missing=0
        for r in rows:
            ch=r.get('hm_chr',r.get('chr_name',''));p=r.get('hm_pos',r.get('chr_position',''))
            if not str(p).isdigit():missing+=1;continue
            state,gt,ids=genome.locus(ch,int(p));ea=r.get('effect_allele','');oa=r.get('other_allele','') or r.get('hm_inferOtherAllele','')
            if state!='called' or len(gt)!=2 or len(ea)!=1 or len(oa)!=1 or set(ea+oa) in (set('AT'),set('CG')) or not set(gt)<=set(ea+oa) or (ch,int(p)) in flagged:missing+=1;continue
            components.append(dict(id=r.get('rsID',''),effect_allele=ea,dosage=gt.count(ea),weight=float(r['effect_weight']),provenance='measured'))
        complete=not missing and bool(rows)
        results.append(dict(id=definition['id'],name=definition.get('name',definition['id']),status='completed' if complete else 'incomplete_coverage',total_variants=len(rows),measured_variants=len(components),unresolved_variants=missing,weighted_sum=sum(r['dosage']*r['weight'] for r in components) if complete else None,source=definition['source'],explanation='A complete weighted effect-allele sum is reported only when every score variant is resolved. It is not a measured trait or an absolute disease-risk estimate.'))
    return dict(status='completed' if definitions else 'reference_not_configured',scores=results)


def region(genome,config,window,work,flagged,targets,seed=2026,masked=False):
    import pysam
    work.mkdir(parents=True,exist_ok=True);ch=str(window['chrom']);calls=[];references={};sample_names=None
    with pysam.VariantFile(window['file']) as vf:
        sample_names=list(vf.header.samples)
        for record in vf:
            if any(len(a)!=1 or a not in 'ACGT' for a in record.alleles):continue
            key=(ch,record.pos)
            if key in targets:
                references[key]=dict(alleles=list(record.alleles),haplotypes={s:list(record.samples[s]['GT']) for s in sample_names})
            if len(record.alleles)!=2 or key in flagged or set(''.join(record.alleles)) in (set('AT'),set('CG')):continue
            state,gt,_=genome.locus(ch,record.pos)
            if state!='called' or len(gt)!=2 or not set(gt)<=set(record.alleles):continue
            hold=int(hashlib.sha256(f'{ch}:{record.pos}'.encode()).hexdigest()[:8],16)%5==0
            calls.append(dict(pos=record.pos,alleles=list(record.alleles),gt=gt,masked=hold))
    if len(calls)<30:return dict(status='insufficient_coverage',input_markers=len(calls),targets={})
    name=f"{window['name']}_{seed}_{'masked' if masked else 'full'}";inp=work/(name+'.vcf');output=work/name
    with inp.open('w') as f:
        f.write(f'##fileformat=VCFv4.2\n##contig=<ID={ch}>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLOCAL_SAMPLE\n')
        for r in calls:
            if masked and r['masked']:continue
            a,b=r['alleles'];g='/'.join(str(int(x==b)) for x in r['gt']);f.write(f'{ch}\t{r["pos"]}\t.\t{a}\t{b}\t.\tPASS\t.\tGT\t{g}\n')
    command=[config['java'],'-Xmx1500m','-jar',config['beagle'],f'gt={inp}',f'ref={window["file"]}',f'map={window["map"]}',f'out={output}','nthreads=1',f'seed={seed}','gp=true','ne=10000']
    log=run_command(command,work,timeout=900);(work/(name+'.stdout.txt')).write_text(log)
    held={r['pos']:r for r in calls if r['masked']};validation=[];out={}
    with pysam.VariantFile(str(output)+'.vcf.gz') as vf:
        for rec in vf:
            sample=rec.samples[0];key=(ch,rec.pos)
            if masked and rec.pos in held and len(rec.alleles)==2:
                truth=held[rec.pos]['gt'].count(rec.alts[0]);gt=sample['GT']
                if None not in gt:validation.append(dict(truth=truth,inferred=sum(gt),gp=float(max(sample.get('GP') or [0]))))
            if key not in targets or key not in references:continue
            gt=sample['GT'];ds=sample.get('DS');ds=list(ds) if isinstance(ds,tuple) else [ds]
            if any(d is None for d in ds) or None in gt:continue
            dosages=[2-sum(ds)]+ds
            out[key]=dict(alleles=list(rec.alleles),dosages=dosages,phased_alleles=[rec.alleles[i] for i in gt],phased=bool(sample.phased),max_gp=float(max(sample.get('GP') or [0])),reference=references[key])
    return dict(status='completed',input_markers=len(calls),targets=out,validation=validation,samples=sample_names)


def imputed_score(genome,config,work,flagged,progress=lambda _:None):
    import numpy as np
    if genome.build!=config.get('build','37'):return dict(status='unsupported_build',explanation='Imputation reference assembly differs from the raw export.')
    score=score_rows(config['score']);targets={(r['hm_chr'],int(r['hm_pos'])):r for r in score};windows=config['windows'];results=[];failures=[]
    def worker(w):return region(genome,config,w,work/w['name'],flagged,targets)
    with ThreadPoolExecutor(max_workers=min(config.get('workers',2),4)) as pool:
        futures={pool.submit(worker,w):w for w in windows}
        try:
            for future in as_completed(futures):
                w=futures[future]
                try:results.append((w,future.result()))
                except Exception as e:failures.append(dict(window=w['name'],reason=str(e)))
                progress(f'Imputation windows {len(results)+len(failures)}/{len(windows)}')
        except BaseException:
            cancel_commands()
            for future in futures:
                future.cancel()
            raise
    calls={};samples=None
    for w,r in sorted(results,key=lambda x:x[0]['name']):
        if r['status']!='completed':continue
        if samples is not None and samples!=r['samples']:raise ValueError('Imputation windows have inconsistent reference samples.')
        samples=r['samples']
        for key,v in r['targets'].items():
            if key in calls:raise ValueError('Overlapping score windows produced duplicate targets.')
            calls[key]=v
    components=[];reference_scores=np.zeros(len(samples or []));unresolved=[]
    for key,row in targets.items():
        v=calls.get(key);ea=row['effect_allele'];oa=row['other_allele'];state,gt,_=genome.locus(*key)
        if not v or ea not in v['alleles'] or oa not in v['alleles']:unresolved.append(dict(id=row['rsID'],reason='reference_or_inference_unavailable'));continue
        direct=state=='called' and len(gt)==2 and set(gt)<={ea,oa} and key not in flagged
        if not direct and v['max_gp']<.99:unresolved.append(dict(id=row['rsID'],reason='imputation_posterior_below_0.99'));continue
        dosage=gt.count(ea) if direct else float(v['dosages'][v['alleles'].index(ea)]);weight=float(row['effect_weight'])
        component=dict(id=row['rsID'],chrom=key[0],pos=key[1],effect_allele=ea,weight=weight,dosage=dosage,provenance='measured' if direct else 'imputed',posterior=None if direct else v['max_gp']);components.append(component)
        ref=v['reference'];idx=ref['alleles'].index(ea)
        for i,sample in enumerate(samples):
            alleles=ref['haplotypes'][sample]
            if None in alleles:reference_scores[i]=np.nan
            else:reference_scores[i]+=weight*alleles.count(idx)
    validation=[]
    for w,r in sorted(results,key=lambda x:x[0]['name'])[:config.get('validation_windows',2)]:
        if r['status']!='completed':continue
        progress('Checking held-out calls in '+w['name'])
        try:validation.extend(region(genome,config,w,work/w['name'],flagged,targets,masked=True)['validation'])
        except Exception as e:failures.append(dict(window=w['name'],reason='Masked validation failed: '+str(e)))
    observed=sum(r['truth']==r['inferred'] for r in validation);complete=len(components)==len(score) and not failures
    total=sum(r['weight']*r['dosage'] for r in components) if complete else None;finite=reference_scores[np.isfinite(reference_scores)]
    return dict(status='completed' if complete else 'incomplete',id=config.get('id','PGS000034'),name=config.get('name','BMI research score'),components=components,unresolved=unresolved,total_variants=len(score),measured=sum(r['provenance']=='measured' for r in components),imputed=sum(r['provenance']=='imputed' for r in components),weighted_sum=total,reference_samples=len(finite),reference_percentile=float(100*np.mean(finite<total)) if complete and len(finite) else None,reference_mean=float(np.mean(finite)) if complete and len(finite) else None,reference_sd=float(np.std(finite,ddof=1)) if complete and len(finite)>1 else None,validation=dict(heldout_calls=len(validation),concordance=observed/len(validation) if validation else None,windows=config.get('validation_windows',2),scope='Deterministically masked measured calls in selected regions; not independent measurements or validation of every imputed site.'),failures=failures,explanation='Published weighted score using separately labelled measured calls and statistically imputed dosages. Its percentile is relative to the supplied reference panel, not a measured BMI, personal disease probability, or clinical threshold.',source=config['source'])


def phase(genome,config,work,flagged):
    if genome.build!=config.get('build','37'):return dict(status='unsupported_build')
    targets={(str(r['chrom']),int(r['pos'])):r for r in config['sites']};runs=[]
    for seed in [2026,17,42]:
        r=region(genome,config,config['window'],work,flagged,targets,seed)
        if r['status']!='completed' or any(k not in r['targets'] for k in targets):return dict(status='insufficient_coverage',explanation='Target phase could not be inferred from this input/reference.')
        patterns=[''.join('1' if r['targets'][k]['phased_alleles'][h]==targets[k]['effect_allele'] else '0' for k in targets) for h in (0,1)]
        runs.append(dict(seed=seed,haplotypes=patterns,posterior={targets[k]['id']:r['targets'][k]['max_gp'] for k in targets},all_phased=all(r['targets'][k]['phased'] for k in targets)))
    direct={v['id']:genome.locus(*k)[1] for k,v in targets.items()};agreement=len({tuple(sorted(r['haplotypes'])) for r in runs})==1
    confident=agreement and all(r['all_phased'] and min(r['posterior'].values())>=.99 for r in runs)
    return dict(status='completed' if confident else 'unresolved',name=config.get('name','Targeted phase'),observed_calls=direct,runs=runs,seed_agreement=agreement,explanation='Each binary string describes the target alleles on one statistically inferred chromosome. Genotype posterior is not phase confidence; agreement across seeds is computational stability. No drug dose or complete metabolizer phenotype is assigned from these sites alone.',source=config['source'])
