import json,collections
from inventory import ROOT,read_raw
rows=read_raw(); bp=collections.defaultdict(list)
for r in rows:bp[(r['chrom'],r['pos'])].append(r)
pcapos=set((p[0],int(p[3])) for l in (ROOT/'references/snpsForPCA.bim').read_text().splitlines() if (p:=l.split()))
panel={p[0]:p[1:3] for l in (ROOT/'references/1000g_panel.tsv').read_text().splitlines()[1:] if (p:=l.split())}
fam=[l.split() for l in (ROOT/'references/omni25/indep.fam').read_text().splitlines()]
keep=[r for r in fam if r[1] in panel]
(ROOT/'analysis/reference_keep.txt').write_text(''.join(f'{r[0]} {r[1]}\n' for r in keep))
candidate=[]; details=[]; excluded=collections.Counter()
for line in (ROOT/'references/omni25/indep.bim').open():
    ch,rs,cm,ps,a1,a2=line.split();p=int(ps)
    if not ch.isdigit() or not 1<=int(ch)<=22:continue
    if (ch,p) not in pcapos:continue
    rr=bp.get((ch,p),[])
    if not rr:excluded['not_in_raw']+=1;continue
    if set(a1+a2) in [set('AT'),set('CG')]:excluded['strand_ambiguous']+=1;continue
    gs={''.join(sorted(r['gt'])) for r in rr if r['gt']!='--'}
    if len(gs)!=1:excluded['missing_or_conflicting']+=1;continue
    gt=next(iter(gs))
    if len(gt)!=2 or not set(gt)<={a1,a2} or a1==a2:excluded['allele_incompatible']+=1;continue
    candidate.append(rs);details.append(dict(id=rs,chrom=ch,pos=p,a1=a1,a2=a2,gt=gt,raw_ids=[r['id'] for r in rr]))
(ROOT/'analysis/ancestry_candidates.txt').write_text('\n'.join(candidate)+'\n')
(ROOT/'analysis/ancestry_marker_mapping.json').write_text(json.dumps(details,indent=2))
print('References total',len(fam),'labelled',len(keep),'populations',collections.Counter(panel[r[1]][0] for r in keep))
print('Candidate markers',len(candidate),'exclusions',excluded)
ys=[]
for (ch,p),rr in sorted(bp.items(),key=lambda x:x[0][1]):
    if ch!='Y':continue
    gs={r['gt'] for r in rr if len(r['gt'])==1 and r['gt'] in 'ACGT'}
    if len(gs)==1:ys.append((p,next(iter(gs))))
(ROOT/'analysis/yhaplo').mkdir(exist_ok=True)
(ROOT/'analysis/yhaplo/sample.genos.txt').write_text('ID\t'+'\t'.join(str(p) for p,g in ys)+'\nSAMPLE\t'+'\t'.join(g for p,g in ys)+'\n')
print('Y unique called nucleotide positions',len(ys))
