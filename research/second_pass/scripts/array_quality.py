"""Compare chip content and published GenomePrep low-quality loci locally."""
from pathlib import Path
import gzip,json,collections,csv
ROOT=Path(__file__).resolve().parents[2];A=ROOT/'second_pass/analysis';R=ROOT/'second_pass/references'
raw=json.loads((ROOT/'analysis/raw_by_id.json').read_text());bp=collections.defaultdict(list)
for r in raw.values():bp[f"{r['chrom']}:{r['pos']}"].append(r)
clusters=collections.defaultdict(set)
for l in gzip.open(R/'chip_clusters.tsv.gz','rt'):
 pos,cc=l.rstrip().split('\t')
 for c in cc.split(','):clusters[c].add(pos)
overlap=[]
for c,sites in clusters.items():
 n=len(sites&bp.keys());overlap.append(dict(cluster=c,panel_sites=len(sites),overlap=n,fraction_of_raw_positions=n/len(bp),fraction_of_cluster=n/len(sites),jaccard=n/len(sites|bp.keys())))
overlap.sort(key=lambda r:-r['jaccard']);best=overlap[0]['cluster']
bad={}
for l in gzip.open(R/'badalleles.tsv.gz','rt'):
 c,sites=l.rstrip().split('\t');bad[c]=set(sites.split(','))
rows=[]
for pos in sorted(bad[best]&bp.keys()):
 for r in bp[pos]:rows.append(dict(**r,quality_flag_cluster=best))
with (A/'array_flagged_calls.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['id','chrom','pos','gt','quality_flag_cluster'],delimiter='\t');w.writeheader();w.writerows(rows)
s=dict(cluster_comparison=overlap,selected_cluster=best,unique_flagged_positions=len(bad[best]&bp.keys()),flagged_raw_rows=len(rows),limits='GenomePrep 2021 population-level low-quality-locus screen, not raw intensity evidence or a clinical assay. Absence from list does not validate a genotype. The best-content cluster is inferred from overlap.')
(A/'array_quality_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
