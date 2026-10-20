"""Local masked-call validation, not a clinical assay or whole-genome imputation."""
from pathlib import Path
import subprocess,json,csv,collections,hashlib,zipfile,os,shutil
ROOT=Path(__file__).resolve().parents[2];A=ROOT/'second_pass/analysis';R=ROOT/'second_pass/references';L=ROOT/'second_pass/logs'
JAVA=os.environ.get('GENETICS_JAVA') or ('/opt/homebrew/opt/openjdk/bin/java' if Path('/opt/homebrew/opt/openjdk/bin/java').exists() else shutil.which('java'))
START,END=30000000,35000000
raw=json.loads((ROOT/'analysis/raw_by_id.json').read_text())
pos=collections.defaultdict(list)
for r in raw.values():
 if r['chrom']=='22' and START<=r['pos']<=END:pos[r['pos']].append(r)
calls={p:rr[0] for p,rr in pos.items() if len({r['gt'] for r in rr})==1 and len(rr[0]['gt'])==2 and set(rr[0]['gt'])<=set('ACGT')}
panel={r['sample']:r['super_pop'] for r in csv.DictReader((ROOT/'references/1000g_panel.tsv').open(),delimiter='\t')}
def prepare():
 counts=collections.Counter();accepted=[]
 err=(L/'unbref.log').open('w')
 proc=subprocess.Popen([JAVA,'-Xmx2g','-jar',str(R/'unbref3.jar'),str(R/'chr22.bref3')],stdout=subprocess.PIPE,stderr=err,text=True)
 for line in proc.stdout:
  if line.startswith('##'):continue
  if line.startswith('#CHROM'):
   samples=line.rstrip().split('\t')[9:];eur=[i for i,s in enumerate(samples) if panel.get(s)=='EUR'];continue
  f=line.rstrip().split('\t',9);p=int(f[1])
  if p<START:continue
  if p>END:break
  counts['reference_region_rows']+=1
  if p not in calls:continue
  r=calls[p];ref,alt=f[3:5]
  if len(ref)!=1 or len(alt)!=1 or not set(r['gt'])<={ref,alt}:counts['nonbiallelic_or_mismatch']+=1;continue
  if set(ref+alt) in [set('AT'),set('CG')]:counts['palindromic_excluded']+=1;continue
  gts=f[9].split('\t');af=sum(gts[i].count('1') for i in eur)/(2*len(eur))
  dos=r['gt'].count(alt);mask=int(hashlib.sha256(str(p).encode()).hexdigest()[:8],16)%5==0
  accepted.append(dict(pos=p,id=r['id'],gt=r['gt'],ref=ref,alt=alt,alt_dosage=dos,eur_af=af,masked=mask))
 proc.terminate();proc.wait();err.close()
 counts.update(raw_region_rows=sum(map(len,pos.values())),usable_raw_positions=len(calls),accepted=len(accepted),masked=sum(r['masked'] for r in accepted),reference_samples=len(samples),eur_reference_samples=len(eur))
 (A/'imputation_input.json').write_text(json.dumps(dict(region=[START,END],counts=counts,variants=accepted),indent=2))
 with (A/'pilot_target.vcf').open('w') as o:
  o.write('##fileformat=VCFv4.2\n##contig=<ID=22,length=51304566>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLOCAL_SAMPLE\n')
  for r in accepted:
   if r['masked']:continue # Entire marker withheld, not a forced missing-ref call.
   g='/'.join(str(int(c==r['alt'])) for c in r['gt'])
   o.write(f"22\t{r['pos']}\t{r['id']}\t{r['ref']}\t{r['alt']}\t.\tPASS\t.\tGT\t{g}\n")
 print(json.dumps(counts,indent=2),flush=True)
def run():
 z=zipfile.ZipFile(R/'plink.GRCh37.map.zip');names=[n for n in z.namelist() if n.endswith('chr22.GRCh37.map')];assert len(names)==1,names
 (R/'chr22.map').write_bytes(z.read(names[0]))
 cmd=[JAVA,'-Xmx4g','-jar',str(R/'beagle.jar'),f'gt={A}/pilot_target.vcf',f'ref={R}/chr22.bref3',f'map={R}/chr22.map',f'out={A}/pilot_imputed',f'chrom=22:{START}-{END}','nthreads=4','seed=2026','gp=true','ne=10000']
 (L/'imputation_command.json').write_text(json.dumps(cmd,indent=2))
 with (L/'imputation_stdout.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
if __name__=='__main__':
 import sys
 if not (A/'imputation_input.json').exists():prepare()
 if '--prepare-only' not in sys.argv:run()
