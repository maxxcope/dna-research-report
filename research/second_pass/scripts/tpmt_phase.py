"""Explore statistical TPMT phase without replacing clinical confirmation."""
from pathlib import Path
import json,collections,subprocess
import pysam
import bmi_imputation as b
ROOT=b.ROOT;A=b.A;R=b.R
w=json.loads((A/'tpmt_reference_manifest.json').read_text())[0]
r=b.worker(w)
sites={18139228:'T',18130918:'C'};phased={};refcalls={}
for rec in pysam.VariantFile(str(ROOT/w['file'])):
 if rec.pos in sites:refcalls[rec.pos]={s:tuple(rec.alleles[a] for a in rec.samples[s]['GT']) for s in rec.samples}
assert len(refcalls)==2
haps=collections.Counter()
for sample in refcalls[18139228]:
 for h in range(2):haps[''.join(str(int(refcalls[p][sample][h]==sites[p])) for p in sites)]+=1
for seed in [2026,17,42]:
 out=b.D/(w['name']+('_full' if seed==2026 else f'_phase_seed{seed}'))
 if seed!=2026:
  cmd=[b.JAVA,'-Xmx2g','-jar',str(R/'beagle.jar'),f'gt={b.D}/{w["name"]}_full.vcf',f'ref={ROOT/w["file"]}',f'map={R}/chr6.map',f'out={out}','nthreads=1',f'seed={seed}','gp=true','ne=10000']
  with (ROOT/'second_pass/logs'/f'tpmt_seed{seed}.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 pc={}
 for rec in pysam.VariantFile(str(out)+'.vcf.gz'):
  if rec.pos in sites:pc[rec.pos]=tuple(rec.alleles[a] for a in rec.samples[0]['GT'])
 pattern=[''.join(str(int(pc[p][h]==sites[p])) for p in sites) for h in range(2)]
 phased[str(seed)]=dict(alleles=pc,variant_haplotypes=pattern,cis=set(pattern)=={'00','11'})
summary=dict(reference_haplotypes=dict(haps),reference_people=len(refcalls[18139228]),phase_runs=phased,input_measured_markers=r['calls'],interpretation='Statistical phasing favors the common cis arrangement if both variants share one inferred chromosome. Repeated seeds and public haplotypes are computational support, not independent assay confirmation; rare trans arrangements and erroneous chip calls remain possible. Clinical TPMT/NUDT15 review before thiopurines remains unchanged.')
(A/'tpmt_phase_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
