"""Download public reference windows around ALL markers in a public BMI score.
This preparation script never reads any personal genotype input.
"""
from pathlib import Path
import gzip,csv,json,collections,concurrent.futures,hashlib,time
import pysam
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'second_pass/references';OUT=R/'bmi_windows';OUT.mkdir(exist_ok=True)
score=list(csv.DictReader((l for l in gzip.open(R/'PGS000034.txt.gz','rt') if not l.startswith('#')),delimiter='\t'))
groups=collections.defaultdict(list)
for r in score:groups[r['hm_chr']].append((max(1,int(r['hm_pos'])-1000000),int(r['hm_pos'])+1000000))
eur=[r['sample'] for r in csv.DictReader((ROOT/'references/1000g_panel.tsv').open(),delimiter='\t') if r['super_pop']=='EUR']
windows={}
for ch,rs in groups.items():
 merged=[]
 for start,end in sorted(rs):
  if merged and start<=merged[-1][1]:merged[-1][1]=max(end,merged[-1][1])
  else:merged.append([start,end])
 windows[ch]=merged
def get(ch):
 url=f'https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr{ch}.1kg.phase3.v5a.vcf.gz'
 v=pysam.VariantFile(url);v.subset_samples(eur);out=[]
 for start,end in windows[ch]:
  name=f'chr{ch}_{start}_{end}';p=OUT/(name+'.vcf.gz');n=0;t=time.time()
  if not p.exists():
   tmp=Path(str(p)+'.part')
   with pysam.VariantFile(str(tmp),'wz',header=v.header) as dest:
    for rec in v.fetch(ch,start-1,end):
     if len(rec.alleles)<2 or any(len(a)!=1 or a not in 'ACGT' for a in rec.alleles):continue
     dest.write(rec);n+=1
   tmp.rename(p)
  out.append(dict(name=name,chrom=ch,start=start,end=end,file=str(p.relative_to(ROOT)),url=url,reference_samples=len(eur),variants_written_this_run=n,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
  print(name,n,round(time.time()-t,1),flush=True)
 return out
if __name__=='__main__':
 out=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  futs={pool.submit(get,ch):ch for ch in windows}
  for fu in concurrent.futures.as_completed(futs):
   try:out.extend(fu.result())
   except Exception as e:print('FAILED',futs[fu],repr(e),flush=True)
   (ROOT/'second_pass/analysis/bmi_reference_manifest.json').write_text(json.dumps(out,indent=2))
 print('completed windows',len(out),'planned',sum(map(len,windows.values())))
