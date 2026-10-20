"""Offline feasibility audit, NOT the application's interpretation engine.

The raw matcher receives no Promethease data. Only its finished results are
compared to the saved report. Output contains aggregates and stays outside the
release directory. Reference snapshots and personal inputs are never uploaded.
"""
import argparse
import base64
from collections import Counter
import hashlib
import html
import json
from pathlib import Path
import re
import sqlite3
import sys
import zlib

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetics_app.parser import parse_genome
from snpedia_rules import evaluate_all


def fields(text):
    """Simple scalar fields for the audit; not a general MediaWiki renderer."""
    return {key.lower(): value.strip() for key, value in
            re.findall(r'^\|\s*([\w]+)\s*=([^\n]*)', text, re.M)}


def normal_text(text):
    return ' '.join(html.unescape(re.sub('<[^>]+>', '', text)).split()).casefold()


def marker_metadata(text):
    match = re.search(r'\{\{\s*(?:Rsnum|23andMe SNP)\b(.*?)\}\}', text, re.S | re.I)
    if not match:
        return {}
    return {key.strip().lower(): value.strip() for key, value in
            re.findall(r'\|\s*([\w]+)\s*=([^|]*)', match[1])}


def report_rows(path):
    chunks = re.findall(
        r"mygenos\.push\.apply\(mygenos,decompressString\('([^']+)'\)\)",
        path.read_text())
    if not chunks:
        raise ValueError('Unsupported saved report format')
    rows = []
    for chunk in chunks:
        rows.extend(json.loads(zlib.decompress(base64.b64decode(chunk))))
    return rows


def oriented_calls(genome, snps, include_chip_aliases=False):
    """Restricted rsID/strand feasibility probe; no position validation.

    Do not use these calls for a production health report. Excludes unresolved
    IDs, conflicting probes, encoded indels, missing orientations and explicit
    build flips. Haploid base calls are doubled only for SNPedia key comparison.
    """
    calls, reasons, counts = {}, {}, Counter()
    proposals = {}
    complement = str.maketrans('ACGT', 'TGCA')
    for key, raw_rows in genome.by_id.items():
        marker = key.lower()
        if marker not in snps:
            reasons[marker] = 'source_marker_absent'
            continue
        metadata = marker_metadata(snps[marker])
        strand = metadata.get('stabilizedorientation', metadata.get('orientation', '')).lower()
        row = raw_rows[0]
        state, genotype, _ = genome.locus(row.chrom, row.pos)
        if key in genome.bad_ids:
            reason = 'conflicting_marker_id'
        elif state != 'called':
            reason = state
        elif metadata.get('flip' + genome.build, '').lower() == 'true':
            reason = 'build_flip_unresolved'
        elif strand not in {'plus', 'minus'}:
            reason = 'orientation_unknown'
        else:
            reason = 'oriented_call'
            genotype = genotype.translate(complement) if strand == 'minus' else genotype
            if len(genotype) == 1:
                genotype *= 2
            calls[marker] = tuple(sorted(genotype))
        reasons[marker] = reason
        counts[reason] += 1
        if include_chip_aliases and re.fullmatch(r'i\d+', marker):
            alias = metadata.get('rsid', '')
            if alias.isdigit() and 'rs' + alias in snps:
                proposals.setdefault('rs' + alias, []).append((key, row))
    for target, probes in proposals.items():
        # A separate experiment: explicit public chip -> rs mappings only.
        # Do not overwrite an existing raw ID or use report redirect fields.
        if target in genome.by_id:
            continue
        metadata = marker_metadata(snps[target])
        strand = metadata.get('stabilizedorientation', metadata.get('orientation', '')).lower()
        observations = set()
        for key, row in probes:
            state, genotype, _ = genome.locus(row.chrom, row.pos)
            if key in genome.bad_ids or state != 'called':
                observations.add(None)
            else:
                observations.add(genotype)
        if len(observations) != 1 or None in observations:
            reasons[target] = 'alias_probe_unresolved'
        elif strand not in {'plus', 'minus'} or metadata.get('flip' + genome.build, '').lower() == 'true':
            reasons[target] = 'alias_orientation_or_build_unresolved'
        else:
            genotype = next(iter(observations))
            genotype = genotype.translate(complement) if strand == 'minus' else genotype
            if len(genotype) == 1:
                genotype *= 2
            calls[target] = tuple(sorted(genotype))
            reasons[target] = 'public_chip_alias_call'
            counts['public_chip_alias_call'] += 1
    return calls, reasons, counts


def audit(raw, database, criteria_path, baseline, include_chip_aliases=False):
    genome = parse_genome(raw.read_bytes())
    with sqlite3.connect(database.resolve().as_uri() + '?mode=ro', uri=True) as connection:
        connection.execute('PRAGMA query_only=ON')
        integrity = connection.execute('PRAGMA quick_check').fetchone()[0]
        snps = dict(connection.execute('SELECT rsid,content FROM snps'))
        genotypes = dict(connection.execute('SELECT id,content FROM genotypes'))
        genosets = dict(connection.execute('SELECT id,content FROM genosets'))
    rules_file = json.loads(criteria_path.read_text())
    if not rules_file.get('complete'):
        raise ValueError('Criteria download is incomplete')
    criteria = {}
    for page in rules_file['pages']:
        if 'revisions' not in page:
            continue
        revision = page['revisions'][0]
        content = revision['slots']['main']['content']
        if hashlib.sha1(content.encode()).hexdigest() != revision['sha1']:
            raise ValueError('Criteria content hash mismatch')
        criteria[page['title'].split('/')[0].lower()] = content

    # Finish all raw inference before loading the comparison report.
    calls, reasons, raw_counts = oriented_calls(genome, snps, include_chip_aliases)
    outcomes, errors = evaluate_all(criteria, calls)
    raw_counts['genotype_page_matches'] = sum(
        marker + '(' + ';'.join(call).lower() + ')' in genotypes
        for marker, call in calls.items())

    rows = report_rows(baseline)
    coverage, agreement, unresolved, unresolved_summaries = (Counter() for _ in range(4))
    baseline_sets = set()
    for row in rows:
        kind = 'variant' if 'rsnum' in row else 'genoset'
        coverage[kind + '_rows'] += 1
        if kind == 'genoset':
            key = row['title'].lower()
            baseline_sets.add(key)
            source = genosets.get(key)
        else:
            marker = row['rsnum'].lower()
            key = marker + row['geno'].lower()
            source = genotypes.get(key)
            coverage['variant_marker_in_snapshot'] += marker in snps
            observed = calls.get(marker)
            if observed is None:
                reason = reasons.get(marker, 'marker_id_absent_from_raw')
                unresolved[reason] += 1
                if row.get('genosummary'):
                    unresolved_summaries[reason] += 1
                agreement['unresolved_or_absent_raw_call'] += 1
            elif '(' + ';'.join(observed).lower() + ')' == row['geno'].lower():
                agreement['exact_report_genotype_agreement'] += 1
                agreement['matching_variant_rows_with_summary'] += bool(row.get('genosummary'))
            else:
                agreement['genotype_disagreement'] += 1
        coverage[kind + '_has_genotype_source'] += source is not None
        source_fields = fields(source or '')
        if row.get('genosummary'):
            coverage[kind + '_with_summary'] += 1
            coverage[kind + '_summary_has_source'] += source is not None
            summary = normal_text(source_fields.get('summary', ''))
            saved = normal_text(row['genosummary'])
            coverage[kind + '_summary_exact_text'] += summary == saved
            coverage[kind + '_summary_equal_after_quote_normalization'] += summary.replace('"', "'") == saved.replace('"', "'")
        if float(row.get('magnitude', 0)) >= 3:
            coverage['magnitude3_rows'] += 1
            coverage['magnitude3_has_source'] += source is not None
    positive = {key for key, value in outcomes.items() if value is True}
    agreement.setdefault('genotype_disagreement', 0)
    return {
        'purpose': 'Feasibility experiment, not production parity certification',
        'public_chip_alias_experiment_enabled': include_chip_aliases,
        'reference': {'sha256': hashlib.sha256(database.read_bytes()).hexdigest(),
                      'sqlite_quick_check': integrity,
                      'snp_pages': len(snps), 'genotypes': len(genotypes), 'genosets': len(genosets)},
        'source_coverage': dict(coverage),
        'raw_matching': dict(raw_counts),
        'baseline_agreement': dict(agreement),
        'unresolved_baseline_rows': dict(unresolved),
        'unresolved_baseline_rows_with_summaries': dict(unresolved_summaries),
        'genoset_experiment': {
            'downloaded_rule_pages': len(criteria),
            'missing_rule_pages': len(rules_file['pages']) - len(criteria),
            'rules_parsed': len(criteria) - len(errors),
            'parse_error_types': dict(Counter(errors.values())),
            'outcomes': dict(Counter(str(value) for value in outcomes.values())),
            'baseline_genosets': len(baseline_sets),
            'baseline_reproduced': len(baseline_sets & positive),
            'baseline_unresolved': sum(outcomes.get(key) is None for key in baseline_sets),
            'baseline_false': sum(outcomes.get(key) is False for key in baseline_sets),
            'extra_positive_genosets': len(positive - baseline_sets)},
        'limits': [
            'One historical report: not proof of present-day Promethease equivalence or medical accuracy.',
            'Source coverage is separate from correct raw genotype matching and explanation rendering.',
            'Raw experiment matches explicit IDs and SNPedia strand metadata, not build-specific coordinates.',
            'Optional chip aliases use explicit public mappings; indels, build flips and missing orientations stay unresolved.',
            'Rule language parses available criteria; three-valued missing-data semantics are our implementation.',
            'Additional genoset matches are unexplained discrepancies, not established new findings.',
            'Genosets with absent rule pages cannot be inferred from their descriptions.',
            'Simple scalar wiki-field extraction is sufficient for this audit, not a full MediaWiki renderer.'
        ]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('raw', 'database', 'criteria', 'promethease', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--include-chip-aliases', action='store_true')
    args = parser.parse_args()
    release = Path(__file__).resolve().parents[1]
    if args.out.resolve().is_relative_to(release):
        parser.error('Write personal comparisons outside the release directory')
    if args.out.resolve() in {p.resolve() for p in (args.raw, args.database, args.criteria, args.promethease)}:
        parser.error('Output must not overwrite an input')
    result = audit(args.raw, args.database, args.criteria, args.promethease, args.include_chip_aliases)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
