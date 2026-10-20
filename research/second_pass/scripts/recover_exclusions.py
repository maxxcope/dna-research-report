"""Audit excluded calls without equating indel codes to DNA sequences."""
from pathlib import Path
import csv,json,collections
ROOT=Path(__file__).resolve().parents[2];A=ROOT/'second_pass/analysis'
raw=json.loads((ROOT/'analysis/raw_by_id.json').read_text())
out=[];counts=collections.Counter()
for r in csv.DictReader((ROOT/'analysis/clinvar_unresolved.tsv').open(),delimiter='\t'):
 gts={''.join(sorted(raw[i]['gt'])) for i in r['raw_ids'].split(';') if set(raw[i]['gt'])<=set('ACGT')}
 reason=r['status']; counts[reason]+=1
 if len(r['ref'])!=1 or len(r['alt'])!=1 or not set(r['ref']+r['alt'])<=set('ACGT'): continue
 if len(gts)!=1:continue
 gt=next(iter(gts))
 if set(gt)<={r['ref'],r['alt']} and reason=='discordant_probes_excluded':
  r.update(recovery='ACGT call retained separately from co-located indel probes',recovered_gt=gt,alt_dosage=gt.count(r['alt']))
  out.append(r)
 elif r['alt'] in gt and not set(gt)<={r['ref'],r['alt']}:
  r.update(recovery='candidate multiallelic ALT presence; other allele requires review',recovered_gt=gt,alt_dosage=gt.count(r['alt']))
  out.append(r)
keys=list(dict.fromkeys(k for r in out for k in r))
with (A/'recovered_clinvar.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(out)
summary=dict(original_exclusions=counts,recovered_records=len(out),recovery_types=collections.Counter(r['recovery'] for r in out),alt_present=sum(r['alt_dosage']>0 for r in out))
(A/'recovery_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
for r in out:
 if r['alt_dosage'] and ('Pathogenic' in r['classification'] or 'Likely_pathogenic' in r['classification']):
  print({k:r[k] for k in ['raw_ids','recovered_gt','ref','alt','gene','classification','variation_id','recovery']})
