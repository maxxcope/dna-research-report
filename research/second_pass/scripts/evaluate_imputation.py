from pathlib import Path
import csv,gzip,json
import numpy as np
ROOT=Path(__file__).resolve().parents[2];A=ROOT/'second_pass/analysis'
truth={r['pos']:r for r in json.loads((A/'imputation_input.json').read_text())['variants'] if r['masked']}
rows=[];total=0
for l in gzip.open(A/'pilot_imputed.vcf.gz','rt'):
 if l.startswith('#'):continue
 f=l.rstrip().split('\t');p=int(f[1]);total+=1
 if p not in truth:continue
 r=truth[p].copy();assert (f[3],f[4])==(r['ref'],r['alt'])
 vals=dict(zip(f[8].split(':'),f[9].split(':')));gp=list(map(float,vals['GP'].split(',')))
 r.update(imputed_dosage=float(vals['DS']),imputed_call=vals['GT'].count('1'),max_posterior=max(gp),prob_true=gp[r['alt_dosage']],correct=vals['GT'].count('1')==r['alt_dosage'])
 rows.append(r)
def describe(rr):
 if not rr:return dict(n=0)
 obs=np.array([r['alt_dosage'] for r in rr]);est=np.array([r['imputed_dosage'] for r in rr]);nonref=[r for r in rr if r['alt_dosage']>0];minor=[r for r in rr if (r['alt_dosage']>0 if r['eur_af']<=.5 else r['alt_dosage']<2)]
 return dict(n=len(rr),concordant=sum(r['correct'] for r in rr),concordance=np.mean([r['correct'] for r in rr]),dosage_r2=float(np.corrcoef(obs,est)[0,1]**2) if len(set(obs))>1 else None,mean_absolute_dosage_error=float(np.mean(abs(obs-est))),reference_only_baseline=float(np.mean(obs==0)),major_homozygote_baseline=float(np.mean([r['alt_dosage']==(0 if r['eur_af']<=.5 else 2) for r in rr])),nonreference_n=len(nonref),nonreference_concordance=float(np.mean([r['correct'] for r in nonref])) if nonref else None,minor_carrier_n=len(minor),minor_carrier_concordance=float(np.mean([r['correct'] for r in minor])) if minor else None)
s=dict(imputed_output_variants=total,withheld_not_returned=sorted(set(truth)-{r['pos'] for r in rows}),all=describe(rows),high_posterior_099=describe([r for r in rows if r['max_posterior']>=.99]),common_eur_maf_ge_005=describe([r for r in rows if min(r['eur_af'],1-r['eur_af'])>=.05]),low_frequency_eur_maf_lt_005=describe([r for r in rows if min(r['eur_af'],1-r['eur_af'])<.05]),limits='One 5 Mb region, one individual; chip-selected held-out calls are comparison truth, not independent clinical truth. Does not validate rare, structural, indel, HLA, star-allele or whole-genome predictions.')
(A/'imputation_validation.json').write_text(json.dumps(s,indent=2))
with (A/'imputation_heldout.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps(s,indent=2))
