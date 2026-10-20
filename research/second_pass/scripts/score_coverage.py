"""Coverage audit of published scores. Missing calls never become reference alleles."""
from pathlib import Path
import json,gzip,csv,collections
ROOT=Path(__file__).resolve().parents[2];A=ROOT/'second_pass/analysis';R=ROOT/'second_pass/references'
raw=json.loads((ROOT/'analysis/raw_by_id.json').read_text());bypos=collections.defaultdict(list)
for r in raw.values():bypos[(r['chrom'],r['pos'])].append(r)
summary=[]
for p in sorted(R.glob('PGS*.txt.gz')):
 counts=collections.Counter();weights=collections.Counter();audit=[]
 meta=json.loads((R/(p.name.split('.')[0]+'.json')).read_text())
 with gzip.open(p,'rt') as f:
  rows=csv.DictReader((l for l in f if not l.startswith('#')),delimiter='\t')
  for r in rows:
   counts['total']+=1;w=abs(float(r['effect_weight']));weights['total_abs_weight']+=w
   ea=r['effect_allele'];oa=r.get('other_allele') or r.get('hm_inferOtherAllele','')
   try:key=(r['hm_chr'],int(r['hm_pos']))
   except (ValueError,KeyError):counts['no_harmonized_position']+=1;continue
   rr=bypos.get(key,[]);status='not_assayed';gt='';ids=''
   if rr:
    ids=';'.join(x['id'] for x in rr);gs={''.join(sorted(x['gt'])) for x in rr if len(x['gt'])==2 and set(x['gt'])<=set('ACGT')}
    if len(gs)==1:
     gt=next(iter(gs))
     if len(ea)!=1 or not set(ea)<=set('ACGT'):status='indel_or_complex_effect_allele'
     elif not set(gt)<=set(ea+oa.replace('/','')):status='allele_mismatch'
     elif len(oa)!=1:status='other_allele_not_uniquely_defined'
     elif set(ea+oa) in [set('AT'),set('CG')]:status='palindromic_orientation_unresolved'
     else:status='usable_direct_call'
    else:status='no_call_or_non_snv'
   counts[status]+=1;weights[status+'_abs_weight']+=w
   if counts['total']<10000 or status!='not_assayed':
    audit.append(dict(score=meta['id'],rsid=r.get('hm_rsID') or r.get('rsID'),chrom=key[0],pos=key[1],effect_allele=ea,other_allele=oa,weight=r['effect_weight'],raw_ids=ids,gt=gt,status=status))
 with (A/(meta['id']+'_coverage.tsv')).open('w') as f:
  wr=csv.DictWriter(f,fieldnames=list(audit[0]),delimiter='\t');wr.writeheader();wr.writerows(audit)
 s=dict(score=meta['id'],trait=meta['trait_reported'],counts=counts,absolute_weight_fraction_direct=weights['usable_direct_call_abs_weight']/weights['total_abs_weight'],direct_fraction=counts['usable_direct_call']/counts['total'],publication=meta['publication'],decision='No full score: incomplete measured coverage; no missing-as-reference substitution or unvalidated proxy score.')
 summary.append(s);print(meta['id'],meta['trait_reported'],dict(counts),s['absolute_weight_fraction_direct'])
(A/'score_coverage_summary.json').write_text(json.dumps(summary,indent=2))
