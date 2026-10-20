"""Local broad GWAS Catalog overlap screen; association is not personal prediction."""
from pathlib import Path
import csv,gzip,json,zipfile,io,re,collections
ROOT=Path(__file__).resolve().parents[2];A=ROOT/'second_pass/analysis';R=ROOT/'second_pass/references'
raw=json.loads((ROOT/'analysis/raw_by_id.json').read_text())
pro_path=ROOT/'analysis/promethease_rows.json'; old_path=ROOT/'analysis/earlier_report.txt'
pro={r.get('rsnum','') for r in json.loads(pro_path.read_text())} if pro_path.exists() else set()
old=old_path.read_text() if old_path.exists() else ''
# Empty comparison inputs do not establish that a marker is absent from Promethease's database.

counts=collections.Counter();stats=collections.defaultdict(lambda:dict(n=0,pmids=set(),traits=set(),best_p=1.0));categories=collections.defaultdict(list)
terms={'sleep':'chronotype|morningness|sleep duration|insomnia','sensory':'taste|bitter|smell|odor|cilantro|hearing|phot.*sneez|asparagus','appearance':'freckl|tanning|hair (color|colour|shape|curvature|morphology)|skin pigment|eye col|baldness','nutrition':'vitamin|carotene|caffeine|lactose|fatty acid|iron (level|status)|ferritin','biomarker':'bilirubin|uric acid|urate|alkaline phosphatase|thyroid|thyrotropin|sex hormone.binding|testosterone|platelet|red blood cell|white blood cell|hemoglobin','immunity':'immunoglobulin|interleukin|C-reactive protein|blood group|coeliac|celiac|microbiome|microbiota','musculoskeletal':'bone mineral|fracture|osteoarthritis|grip strength|creatine kinase','reproduction':'menopause|menarche|fertility|sperm|azoospermia'}
regs={k:re.compile(v,re.I) for k,v in terms.items()}
z=zipfile.ZipFile(R/'gwas_catalog.zip')
fields=['raw_id','raw_gt','raw_chrom','raw_pos37','PUBMEDID','DATE','DISEASE/TRAIT','MAPPED_GENE','STRONGEST SNP-RISK ALLELE','P-VALUE','OR or BETA','95% CI (TEXT)','INITIAL SAMPLE SIZE','REPLICATION SAMPLE SIZE','STUDY ACCESSION','CONTEXT','in_promethease','in_earlier_pdf']
with gzip.open(A/'gwas_overlaps.tsv.gz','wt') as o:
 w=csv.DictWriter(o,fieldnames=fields,delimiter='\t');w.writeheader()
 for r in csv.DictReader(io.TextIOWrapper(z.open(z.namelist()[0])),delimiter='\t'):
  counts['catalog_rows']+=1
  ids=set(re.findall(r'rs\d+',r['SNPS']+' '+r['STRONGEST SNP-RISK ALLELE']))
  if r['SNP_ID_CURRENT'].isdigit():ids.add('rs'+r['SNP_ID_CURRENT'])
  hits=ids&raw.keys()
  if not hits:continue
  counts['catalog_rows_with_raw_id_overlap']+=1
  try:p=float(r['P-VALUE'])
  except ValueError:continue
  for rs in sorted(hits):
   rr=raw[rs]
   d={k:r.get(k,'') for k in fields};d.update(raw_id=rs,raw_gt=rr['gt'],raw_chrom=rr['chrom'],raw_pos37=rr['pos'],in_promethease=rs in pro,in_earlier_pdf=rs in old)
   w.writerow(d)
   if rr['gt']=='--' or p>5e-8:continue
   counts['significant_called_overlap_records']+=1
   s=stats[rs];s['n']+=1;s['pmids'].add(r['PUBMEDID']);s['traits'].add(r['DISEASE/TRAIT']);s['best_p']=min(s['best_p'],p)
   for cat,rx in regs.items():
    if rx.search(r['DISEASE/TRAIT']):categories[cat].append(d)
out=[]
for rs,s in stats.items():
 out.append(dict(rs=rs,gt=raw[rs]['gt'],records=s['n'],publications=len(s['pmids']),traits=len(s['traits']),best_p=s['best_p'],in_promethease=rs in pro,in_earlier_pdf=rs in old))
out.sort(key=lambda r:(-r['publications'],r['best_p']))
with (A/'gwas_marker_summary.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
for cat,rr in categories.items():
 with gzip.open(A/(cat+'_gwas_candidates.tsv.gz'),'wt') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rr)
counts['distinct_significant_called_markers']=len(stats);counts['distinct_significant_markers_absent_promethease']=sum(not r['in_promethease'] for r in out)
(A/'gwas_scan_summary.json').write_text(json.dumps(dict(counts=counts,category_overlap_records={k:len(v) for k,v in categories.items()},limits='rsID overlap is discovery only. Risk/effect allele orientation, LD, phenotype units and publication evidence must be reviewed before interpreting a finding. Counts are not independent effects or traits predicted.'),indent=2))
print(json.dumps(counts,indent=2))
for cat,rr in categories.items():
 dd=collections.defaultdict(list)
 for r in rr:dd[r['raw_id']].append(r)
 ss=sorted(dd,key=lambda rs:(-len(set(x['PUBMEDID'] for x in dd[rs])),min(float(x['P-VALUE']) for x in dd[rs])))[:12]
 print('\n',cat,len(rr))
 for rs in ss:
  r=min(dd[rs],key=lambda x:float(x['P-VALUE']));print(rs,r['raw_gt'],r['MAPPED_GENE'],r['DISEASE/TRAIT'],r['P-VALUE'],r['STRONGEST SNP-RISK ALLELE'],'PMID',r['PUBMEDID'])
