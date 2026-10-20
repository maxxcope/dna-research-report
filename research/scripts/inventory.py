from pathlib import Path
import json, csv, re, zlib, base64, hashlib, os
from collections import Counter, defaultdict

ROOT=Path(__file__).resolve().parents[1]
# Explicit local inputs only; never fall back to a personal workspace.
RAW=Path(os.environ['GENETICS_RAW']) if os.environ.get('GENETICS_RAW') else ROOT/'examples/synthetic_genotypes.txt'
PROM=Path(os.environ['GENETICS_PROMETHEASE']) if os.environ.get('GENETICS_PROMETHEASE') else None
PDF=Path(os.environ['GENETICS_PRIOR_PDF']) if os.environ.get('GENETICS_PRIOR_PDF') else None
def read_raw():
    rows=[]
    for line in RAW.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'): continue
        rs,ch,p,gt=line.split()
        rows.append(dict(id=rs,chrom=ch,pos=int(p),gt=gt))
    return rows

def main():
    (ROOT/'analysis').mkdir(exist_ok=True)
    rows=read_raw(); byid={r['id']:r for r in rows}
    bypos=defaultdict(list)
    counts=defaultdict(Counter)
    for r in rows:
        bypos[(r['chrom'],r['pos'])].append(r)
        c=counts[r['chrom']]; c['rows']+=1
        c['missing']+=r['gt']=='--'
        c['called']+=r['gt']!='--'
        c['heterozygous_acgt']+=len(r['gt'])==2 and set(r['gt'])<=set('ACGT') and r['gt'][0]!=r['gt'][1]
        c['acgt_diploid']+=len(r['gt'])==2 and set(r['gt'])<=set('ACGT')
        c['indel_coded']+=bool(set(r['gt'])&set('ID'))
    dup=[v for v in bypos.values() if len(v)>1]
    conflict=[v for v in dup if len({''.join(sorted(r['gt'])) for r in v if r['gt']!='--'})>1]
    qc=dict(total_rows=len(rows),unique_ids=len(byid),chromosomes=counts,duplicate_positions=len(dup),conflicting_positions=conflict,
        missing=sum(r['gt']=='--' for r in rows),called=sum(r['gt']!='--' for r in rows),
        genotype_counts=Counter(r['gt'] for r in rows),
        sources=[dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [RAW,PROM,PDF] if p is not None])
    (ROOT/'analysis/qc.json').write_text(json.dumps(qc,indent=2))
    with (ROOT/'analysis/raw_by_id.json').open('w') as f:json.dump(byid,f)
    s=PROM.read_text() if PROM is not None else ''
    chunks=re.findall(r"mygenos\.push\.apply\(mygenos,decompressString\('([^']+)'\)\)",s)
    pro=[]
    for chunk in chunks: pro.extend(json.loads(zlib.decompress(base64.b64decode(chunk))))
    (ROOT/'analysis/promethease_rows.json').write_text(json.dumps(pro,indent=2))
    print(json.dumps({k:v for k,v in qc.items() if k not in ['sources','conflicting_positions']},indent=2))
    print('Promethease:',len(chunks),'chunks;',len(pro),'rows')
    print('Duplicate conflicts:',len(conflict))

if __name__=='__main__':main()
