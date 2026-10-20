"""Local GRCh37 exact-position, exact-allele SNV comparison; no imputation."""
from pathlib import Path
import csv,gzip,json
from collections import defaultdict,Counter
from inventory import read_raw,ROOT

rows=read_raw(); bypos=defaultdict(list)
for r in rows:bypos[(r['chrom'],r['pos'])].append(r)
stats=Counter(); hits=[]; overlaps=[]; headers=[]; unresolved=[]
for line in gzip.open(ROOT/'references/clinvar_GRCh37_2026-10-12.vcf.gz','rt'):
    if line.startswith('#'):
        headers.append(line);continue
    stats['vcf_records']+=1
    ch,pos,vid,ref,alt,qual,filt,info=line.rstrip().split('\t')[:8]
    matched=bypos.get((ch,int(pos)))
    if not matched:continue
    stats['position_overlap_records']+=1
    inf=dict(x.split('=',1) for x in info.split(';') if '=' in x)
    d=dict(chrom=ch,pos=pos,variation_id=vid,ref=ref,alt=alt,raw_ids=';'.join(r['id'] for r in matched),
        genotypes=';'.join(r['gt'] for r in matched),classification=inf.get('CLNSIG',''),review=inf.get('CLNREVSTAT',''),
        gene=inf.get('GENEINFO',''),condition=inf.get('CLNDN',''),rs=inf.get('RS',''),
        af_tgp=inf.get('AF_TGP',''),af_exac=inf.get('AF_EXAC',''),origin=inf.get('ORIGIN',''),
        conflict=inf.get('CLNSIGCONF',''),hgvs=inf.get('CLNHGVS',''))
    gts={''.join(sorted(r['gt'])) for r in matched if r['gt']!='--'}
    if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':
        d['status']='non_SNV_not_inferred';unresolved.append(d);stats['non_snv_overlap']+=1;continue
    if len(gts)>1:
        d['status']='discordant_probes_excluded';unresolved.append(d);stats['discordant_snv_excluded']+=1;continue
    if not gts:
        d['status']='no_call';stats['no_call_records']+=1;overlaps.append(d);continue
    gt=next(iter(gts))
    if not set(gt)<={ref,alt}:
        d['status']='allele_mismatch_excluded';unresolved.append(d);stats['allele_mismatch_excluded']+=1;continue
    d['status']='ALT_present' if alt in gt else 'REF_only'
    d['alt_dosage']=gt.count(alt)
    overlaps.append(d);stats[d['status']]+=1
    if alt in gt:hits.append(d)

def write(name,data):
    keys=list(dict.fromkeys(k for r in data for k in r))
    with (ROOT/'analysis'/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(data)
write('clinvar_snv_overlaps.tsv',overlaps);write('clinvar_alt_matches.tsv',hits);write('clinvar_unresolved.tsv',unresolved)
flags=[r for r in hits if any(t in r['classification'].lower() for t in ['pathogenic','risk','drug','protective','affects'])]
write('clinvar_review_queue.tsv',flags)
(ROOT/'analysis/clinvar_headers.txt').write_text(''.join(headers))
stats['alt_classifications']=Counter(r['classification'] for r in hits)
(ROOT/'analysis/clinvar_summary.json').write_text(json.dumps(stats,indent=2))
print(json.dumps(stats,indent=2))
for r in flags:
    if 'pathogenic' in r['classification'].lower() and 'conflicting' not in r['classification'].lower():
        print({k:r[k] for k in ['raw_ids','genotypes','gene','classification','review','condition','variation_id']})
