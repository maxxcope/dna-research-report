"""Complete a small public BMI score locally; retain measured vs imputed provenance.
This is a research score, never a BMI measurement or calibrated clinical risk.
"""
from pathlib import Path
import json,csv,gzip,collections,subprocess,concurrent.futures,hashlib,zipfile,os,shutil
import numpy as np,pysam
ROOT=Path(__file__).resolve().parents[2];A=ROOT/'second_pass/analysis';R=ROOT/'second_pass/references';L=ROOT/'second_pass/logs';D=A/'bmi';D.mkdir(exist_ok=True)
JAVA=os.environ.get('GENETICS_JAVA') or ('/opt/homebrew/opt/openjdk/bin/java' if Path('/opt/homebrew/opt/openjdk/bin/java').exists() else shutil.which('java'))
raw=json.loads((ROOT/'analysis/raw_by_id.json').read_text());bypos=collections.defaultdict(list)
for r in raw.values():bypos[(r['chrom'],r['pos'])].append(r)
bad={(r['chrom'],int(r['pos'])) for r in csv.DictReader((A/'array_flagged_calls.tsv').open(),delimiter='\t')}
score=list(csv.DictReader((l for l in gzip.open(R/'PGS000034.txt.gz','rt') if not l.startswith('#')),delimiter='\t'));scorekeys={(r['hm_chr'],int(r['hm_pos'])):r for r in score}
windows=json.loads((A/'bmi_reference_manifest.json').read_text());z=zipfile.ZipFile(R/'plink.GRCh37.map.zip')
for ch in {w['chrom'] for w in windows}:
 n=[n for n in z.namelist() if n.endswith(f'chr{ch}.GRCh37.map')];assert len(n)==1
 (R/f'chr{ch}.map').write_bytes(z.read(n[0]))
def worker(w):
 ch=w['chrom'];name=w['name'];ref=ROOT/w['file'];calls=[];score_refs={}
 vf=pysam.VariantFile(str(ref));samples=list(vf.header.samples)
 for rec in vf:
  p=rec.pos;key=(ch,p);rr=bypos.get(key,[])
  if key in scorekeys:
   sr=scorekeys[key];ea=sr['effect_allele'];oa=sr['other_allele']
   if {ea,oa}<=set(rec.alleles):score_refs[p]=dict(alleles=rec.alleles,reference_effect_dosages=[sum(rec.alleles[a]==ea for a in rec.samples[s]['GT']) for s in samples])
  if len(rec.alleles)!=2:continue
  refa,alt=rec.alleles
  if key in bad or set(refa+alt) in [set('AT'),set('CG')]:continue
  gts={''.join(sorted(r['gt'])) for r in rr if len(r['gt'])==2 and set(r['gt'])<=set('ACGT')}
  if len(gts)!=1:continue
  gt=next(iter(gts))
  if not set(gt)<={refa,alt}:continue
  af=sum(sum(rec.samples[s]['GT']) for s in samples)/(2*len(samples))
  calls.append(dict(pos=p,gt=gt,ref=refa,alt=alt,dosage=gt.count(alt),af=af,masked=int(hashlib.sha256((ch+':'+str(p)).encode()).hexdigest()[:8],16)%5==0))
 for mode in ['full','masked']:
  inp=D/(name+'_'+mode+'.vcf');out=D/(name+'_'+mode)
  with inp.open('w') as f:
   f.write(f'##fileformat=VCFv4.2\n##contig=<ID={ch}>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLOCAL_SAMPLE\n')
   for r in calls:
    if mode=='masked' and r['masked']:continue
    gt='/'.join(str(int(a==r['alt'])) for a in r['gt'])
    f.write(f"{ch}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t.\tGT\t{gt}\n")
  if not Path(str(out)+'.vcf.gz').exists():
   cmd=[JAVA,'-Xmx2g','-jar',str(R/'beagle.jar'),f'gt={inp}',f'ref={ref}',f'map={R}/chr{ch}.map',f'out={out}','nthreads=1','seed=2026','gp=true','ne=10000']
   with (L/(name+'_'+mode+'.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 held={r['pos']:r for r in calls if r['masked']};validation=[]
 for rec in pysam.VariantFile(str(D/(name+'_masked.vcf.gz'))):
  if rec.pos not in held:continue
  r=held[rec.pos];s=rec.samples[0];ds=s['DS'];ds=ds[0] if isinstance(ds,tuple) else ds
  validation.append(dict(chrom=ch,pos=rec.pos,truth=r['dosage'],imputed=sum(s['GT']),dosage=float(ds),gp=max(s['GP']),eur_af=r['af']))
 srout=[]
 for rec in pysam.VariantFile(str(D/(name+'_full.vcf.gz'))):
  p=rec.pos;key=(ch,p)
  if key not in scorekeys or p not in score_refs:continue
  sr=scorekeys[key];s=rec.samples[0];ds=s['DS'];ds=list(ds) if isinstance(ds,tuple) else [ds]
  ea_index=rec.alleles.index(sr['effect_allele']);ds=float(ds[ea_index-1]) if ea_index else 2-sum(ds)
  rr=bypos.get(key,[]);direct=[r for r in rr if len(r['gt'])==2 and set(r['gt'])<={sr['effect_allele'],sr['other_allele']}]
  # Historical scoring logic; palindromic orientation requires explicit validation
  # before this method can be generalized to other samples.
  use_direct=bool(direct) and key not in bad
  srout.append(dict(rsid=sr['rsID'],chrom=ch,pos=p,effect_allele=sr['effect_allele'],other_allele=sr['other_allele'],weight=float(sr['effect_weight']),dosage=direct[0]['gt'].count(sr['effect_allele']) if use_direct else ds,provenance='measured' if use_direct else 'imputed',raw_gt=direct[0]['gt'] if direct else '',imputed_effect_dosage=ds,imputation_max_gp=float(max(s['GP'])),reference_effect_dosages=score_refs[p]['reference_effect_dosages']))
 print(name,'inputs',len(calls),'validation',len(validation),'score_sites',len(srout),flush=True)
 return dict(window=name,calls=len(calls),validation=validation,score=srout,samples=samples)
if __name__=='__main__':
 results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for result in pool.map(worker,windows):results.append(result)
 (A/'bmi_imputation_results.json').write_text(json.dumps(results))
 ss=[r for x in results for r in x['score']];ss.sort(key=lambda r:(int(r['chrom']),r['pos']));assert len({r['rsid'] for r in ss})==len(ss), 'duplicate score variants'
 val=[r for x in results for r in x['validation']]
 for name,rows in [('bmi_score_components.tsv',ss),('bmi_masked_validation.tsv',val)]:
  fields=[k for k in rows[0] if k!='reference_effect_dosages']
  with (A/name).open('w') as f:
   wr=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');wr.writeheader();wr.writerows(rows)
 def describe(rows):
  o=np.array([r['truth'] for r in rows]);d=np.array([r['dosage'] for r in rows]);called=np.array([r['imputed'] for r in rows]);minor=np.array([r['truth']>0 if r['eur_af']<=.5 else r['truth']<2 for r in rows])
  return dict(n=len(rows),concordance=float(np.mean(o==called)),dosage_r2=float(np.corrcoef(o,d)[0,1]**2),minor_carriers=int(minor.sum()),minor_carrier_concordance=float(np.mean(o[minor]==called[minor])),major_homozygote_baseline=float(np.mean([r['truth']==(0 if r['eur_af']<=.5 else 2) for r in rows])))
 result=dict(complete=len(ss)==97,components=len(ss),measured=sum(r['provenance']=='measured' for r in ss),imputed=sum(r['provenance']=='imputed' for r in ss),low_posterior_imputed=[r['rsid'] for r in ss if r['provenance']=='imputed' and r['imputation_max_gp']<.99],validation=describe(val),validation_posterior_ge099=describe([r for r in val if r['gp']>=.99]))
 if len(ss)==97:
  us=sum(r['weight']*r['dosage'] for r in ss);refs=np.sum([r['weight']*np.array(r['reference_effect_dosages']) for r in ss],axis=0);samples=results[0]['samples'];assert all(x['samples']==samples for x in results)
  panel={r['sample']:r['pop'] for r in csv.DictReader((ROOT/'references/1000g_panel.tsv').open(),delimiter='\t')};labs=np.array([panel[s] for s in samples]);north=np.isin(labs,['CEU','GBR'])
  result.update(weighted_sum=us,publication_normalized_score=us*97/sum(r['weight'] for r in ss),eur_mean=float(refs.mean()),eur_sd=float(refs.std(ddof=1)),eur_empirical_percentile=float(100*np.mean(refs<us)),eur_n=len(refs),ceu_gbr_percentile=float(100*np.mean(refs[north]<us)),ceu_gbr_n=int(north.sum()),score_definition='Sum of PGS000034 effect weight times effect-allele dosage; publication normalization multiplies by 97/sum(weights). Empirical rank in a small public genetic reference, not a measured BMI or a calibrated clinical prediction.')
  pdrows=[dict(sample=s,population=panel[s],score=float(v)) for s,v in zip(samples,refs)]
  with (A/'bmi_reference_scores.tsv').open('w') as f:
   wr=csv.DictWriter(f,fieldnames=list(pdrows[0]),delimiter='\t');wr.writeheader();wr.writerows(pdrows)
 (A/'bmi_score_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
