from inventory import ROOT,read_raw
from collections import defaultdict
import json
refpath=ROOT/'references/haplogrep3/data/phylotrees/phylotree-rcrs/17.3/rcrs.fasta'
if not refpath.exists():
    refpath=next((ROOT/'references/haplogrep3').rglob('rcrs.fasta'))
ref=''.join(l.strip() for l in refpath.read_text().splitlines() if not l.startswith('>'))
bp=defaultdict(list)
for r in read_raw():
    if r['chrom']=='MT':bp[r['pos']].append(r)
called={};excluded={}
for p,rr in bp.items():
    gs={r['gt'] for r in rr if r['gt']!='--'}
    if len(gs)==1 and next(iter(gs)) in ['A','C','G','T'] and p<=len(ref):called[p]=next(iter(gs))
    else:excluded[p]=list(gs)
diff=[f'{p}{a}' for p,a in sorted(called.items()) if a!=ref[p-1]]
ranges=';'.join(str(p) for p in sorted(called))
(ROOT/'analysis/mito.hsd').write_text('SAMPLE\t"'+ranges+'"\t\t'+'\t'.join(diff)+'\n')
(ROOT/'analysis/mito_input_summary.json').write_text(json.dumps(dict(reference=str(refpath),reference_length=len(ref),called_unique=len(called),missing_or_ambiguous=excluded,differences=diff),indent=2))
print('MT',len(called),'unique callable positions;',len(diff),'rCRS differences:',diff)
