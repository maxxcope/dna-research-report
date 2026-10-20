"""Optional LOCAL assessment only. Never used by the application or its reports.

Writes aggregate coverage counts, not genotype rows or copied SNPedia prose.
Use an output path outside the release folder for any personal comparison.
"""
import argparse
import base64
from collections import Counter
import json
from pathlib import Path
import re
import zlib

def interpretation_comparison(rows,report):
    knowledge=report.get('knowledge',{})
    if not knowledge:return {'status':'broad_layer_not_present'}
    findings={r['id']:r for r in knowledge['findings']};unresolved={r['id']:r['reason'] for r in knowledge['unresolved']}
    aliases={alias:r['id'] for r in findings.values() for alias in r.get('public_aliases',[])}
    counts=Counter();reasons=Counter();baseline_sets=set()
    for row in rows:
        if 'rsnum' in row:
            key=row['rsnum'].lower();result=findings.get(key);counts['baseline_variant_rows']+=1
            if not result:
                if key in aliases:
                    counts['represented_by_explicit_public_alias']+=1
                    continue  # Do not assume the saved chip page uses the canonical page's strand notation.
                counts['unresolved_baseline_variants']+=1;reasons[unresolved.get(key,'not_resolved_from_raw_input')]+=1;continue
            if result.get('source_genotype','').lower()!=row['geno'].lower():counts['genotype_disagreements']+=1;continue
            counts['exact_genotype_matches']+=1
            if row.get('genosummary'):counts['baseline_summaries_with_exact_genotype_match']+=1
        else:baseline_sets.add(row['title'].lower())
    actual_sets={r['id'] for r in knowledge['findings'] if r['kind']=='genoset'}
    rule_status={r['id']:r['status'] for r in knowledge['genoset_audit']}
    counts['baseline_genosets']=len(baseline_sets);counts['matched_baseline_genosets']=len(baseline_sets&actual_sets)
    counts['extra_genosets']=len(actual_sets-baseline_sets);counts.setdefault('genotype_disagreements',0)
    return dict(counts=dict(counts),unresolved_reasons=dict(reasons),unmatched_genoset_statuses=dict(Counter(rule_status.get(k,'absent') for k in baseline_sets-actual_sets)),scope='Comparison of independently generated results with a historical saved report. Extra findings are differences, not evidence of superiority.')

def compare(html,report):
    chunks=re.findall(r"mygenos\.push\.apply\(mygenos,decompressString\('([^']+)'\)\)",html)
    if not chunks:raise ValueError('Unsupported Promethease HTML format: no compressed records found.')
    rows=[]
    for c in chunks:rows.extend(json.loads(zlib.decompress(base64.b64decode(c))))
    baseline={str(r.get('rsnum','')).lower() for r in rows if re.fullmatch(r'rs\d+',str(r.get('rsnum','')).lower())}
    cv={r['rsid'] for r in report['clinvar'] if r['rsid']}
    traits={r['rsid'] for r in report['traits'] if r['status']=='called'}
    gwas={r['rs'] for r in report['gwas']}
    return dict(
        generated_findings=interpretation_comparison(rows,report),
        baseline=dict(entries=len(rows),distinct_rsids=len(baseline),genoset_entries=sum(str(r.get('title','')).startswith('gs') for r in rows),
                      nonempty_genotype_summaries=sum(bool(r.get('genosummary')) for r in rows),with_clinvar_field=sum(bool(r.get('clinvar_1')) for r in rows),
                      with_population_frequency=sum(bool(r.get('popfreqallgenos')) for r in rows),with_reference_count=sum(bool(r.get('numrefs')) for r in rows),with_topics=sum(bool(r.get('topic')) for r in rows)),
        app=dict(clinvar_annotation_records=len(report['clinvar']),clinvar_distinct_rsids=len(cv),interpreted_curated_traits=len(traits),gwas_overlap_records=report['gwas_counts'].get('total_overlap_records',0),gwas_retained_distinct_rsids=len(gwas)),
        overlap=dict(baseline_rsids_with_clinvar_alt_annotation=len(baseline&cv),baseline_rsids_with_curated_interpretation=len(baseline&traits),baseline_rsids_with_retained_gwas_lead=len(baseline&gwas),clinvar_rsids_absent_from_saved_baseline=len(cv-baseline)),
        limits=[
            'One saved report and one raw export; not a benchmark of current Promethease or clinical accuracy.',
            'Overlap compares marker IDs, not genotype-specific interpretation agreement. Counts are not directly interchangeable.',
            'ClinVar annotations include benign and uncertain records and are not all clinically useful findings.',
            'GWAS leads are capped in the app and are not genotype interpretations.',
            'No independent medical validation or assertion of equivalence/superiority is made.'
        ])

def main():
    p=argparse.ArgumentParser();p.add_argument('--promethease',type=Path,required=True);p.add_argument('--report-json',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    repo=Path(__file__).resolve().parents[1]
    if a.out.resolve().is_relative_to(repo):p.error('Write personal comparisons outside the release folder.')
    data=compare(a.promethease.read_text(),json.loads(a.report_json.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
if __name__=='__main__':main()
