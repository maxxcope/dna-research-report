import json
from pathlib import Path
import re
import tempfile
import unittest

from genetics_app.bundle import at_path, compact, expand_part
from genetics_app.retrieval import EvidenceIndex
from genetics_app.report import export_report
from tests.test_bundle import example, part_data, rebuild
from tests.test_navigation import links


def fixture():
    r = example()
    r['knowledge'] = {'findings': [
        {'id': 'gs101', 'kind': 'genoset', 'summary': 'EXAMPLE source pattern', 'dependencies': ['rs101'],
         'assay_warnings': ['Keep the source caution']},
        {'id': 'rs101', 'kind': 'variant', 'gene': 'EXAMPLE', 'raw_ids': ['i101'], 'summary': 'Measured support'},
        {'id': 'gs999', 'kind': 'genoset', 'summary': 'Different subject', 'dependencies': ['rs101', 'rs999']},
        {'id': 'rs999', 'kind': 'variant', 'gene': 'UNRELATED', 'summary': 'Never follow the reverse rule'},
        {'id': 'rs1010', 'kind': 'variant', 'gene': 'EXAMPLE2', 'summary': 'Not an exact search match'}],
        'unresolved': [], 'genoset_audit': [
            {'id': 'gs101', 'rule': 'rs101(A;A)', 'status': 'matched'},
            {'id': 'gs999', 'rule': 'and(rs101(A;A),rs999(G;G))', 'status': 'not_matched'}]}
    r['finding_audit'] = {'flagged_findings': [{'id': 'i101', 'flags': ['Probe warning']}],
                          'unresolved_variants': [], 'genosets': r['knowledge']['genoset_audit']}
    r['clinvar'] = [{'rsid': 'rs101', 'classification': 'conflicting', 'conditions': 'Example condition'}]
    r['conclusions'] = {'cards': [{'id': 'medication:EXAMPLE', 'title': 'EXAMPLE source patterns',
                                  'bucket': 'medication', 'source_patterns': [{'id': 'gs101'}],
                                  'limits': 'Not a complete medication phenotype. No dose recommendation.'}],
                        'evidence_limits': 'No independent scientific validation inferred.'}
    return r


class RetrievalTests(unittest.TestCase):
    def test_forward_support_and_probe_warning_without_reverse_rule_or_substring(self):
        r = fixture();index = EvidenceIndex(r)
        hits = index.search('EXAMPLE')
        self.assertEqual(hits[0]['path'], '/conclusions/cards/0')
        self.assertNotIn('/knowledge/findings/4', [h['path'] for h in hits])
        result = index.select([hits[0]['path']])
        self.assertEqual(set(result['records']), {
            '/conclusions/cards/0', '/knowledge/findings/0', '/knowledge/findings/1',
            '/knowledge/genoset_audit/0', '/finding_audit/flagged_findings/0', '/clinvar/0',
            '/limits', '/conclusions/evidence_limits'})
        for path, value in result['records'].items():
            self.assertEqual(value, at_path(r, path))
        self.assertNotIn('UNRELATED', compact(result))

    def test_focused_files_have_finite_closed_plan_and_only_expected_records(self):
        r = fixture()
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(r, tmp)['bundle']
            original, manifest, _ = rebuild(folder)
            self.assertEqual(original, r)
            view = manifest['focus']['/conclusions/cards/0']
            rows = [row for file in view['files'] for row in expand_part(part_data(folder / file))]
            self.assertEqual({row['path'] for row in rows}, set(EvidenceIndex(r).select(['/conclusions/cards/0'])['records']))
            closed = links(folder, view['map'])
            self.assertEqual(closed, {(folder / f).resolve() for f in [view['map'], *view['files']]})
            for file in view['files']:
                self.assertIn('Then stop', (folder / file).read_text())
                self.assertNotIn('Use index.tsv', (folder / file).read_text())
                self.assertNotIn('UNRELATED', (folder / file).read_text())
                self.assertLessEqual((folder / file).stat().st_size, manifest['max_part_bytes'])
            self.assertIn((folder / view['map']).resolve(), links(folder, 'index.md'))

    def test_missing_forward_record_stated_without_claiming_missing_dna(self):
        r = fixture();r['knowledge']['findings'][0]['dependencies'].append('rs404')
        result = EvidenceIndex(r).select(['/conclusions/cards/0'])
        self.assertIn('rs404', result['no_separate_record_for'])
        self.assertIn('not absence of a genotype', result['stop'])

    def test_conditional_and_assay_summaries_keep_their_named_source_evidence(self):
        r = fixture()
        r['conclusions']['cards'] += [
            {'id': 'conditional:gs999', 'title': 'Conditional example', 'bucket': 'conditional',
             'observed_support': {'rs101': {'genotype': 'AA'}}, 'limits': 'Incomplete rule'},
            {'id': 'assay:rs101', 'title': 'Assay caution', 'bucket': 'assay_caution'}]
        index = EvidenceIndex(r)
        conditional = index.select(['/conclusions/cards/1'])['records']
        self.assertIn('/knowledge/genoset_audit/1', conditional)
        self.assertIn('/knowledge/findings/3', conditional)
        assay = index.select(['/conclusions/cards/2'])['records']
        self.assertIn('/knowledge/findings/1', assay)
        self.assertIn('/finding_audit/flagged_findings/0', assay)
        self.assertNotIn('/knowledge/genoset_audit/1', assay)

    def test_invalid_scope_and_empty_query_fail_clearly(self):
        index = EvidenceIndex(fixture())
        for paths in ([], ['/knowledge/findings/999999'], ['/limits']):
            with self.assertRaises(ValueError):index.select(paths)
        with self.assertRaises(ValueError):index.search('  ')
        self.assertEqual(index.search('NO_SUCH_GENE'), [])

    def test_medication_summary_keeps_other_gene_alleles_and_unresolved_calls(self):
        r = fixture()
        r['knowledge']['findings'].append({'id': 'rs202', 'gene': 'EXAMPLE',
                                          'summary': 'Different allele has a different source label'})
        r['knowledge']['unresolved'] = [{'id': 'rs203', 'gene': 'EXAMPLE', 'reason': 'Missing allele definition'}]
        r['clinvar'].append({'rsid': 'rs204', 'gene': 'EXAMPLE:123|SECOND:456', 'classification': 'other'})
        r['finding_audit']['flagged_findings'].append({'id': 'rs202', 'flags': ['Keep this caution']})
        selected = EvidenceIndex(r).select(['/conclusions/cards/0'])['records']
        for path in ('/knowledge/findings/5', '/knowledge/unresolved/0', '/clinvar/1', '/finding_audit/flagged_findings/1'):
            self.assertIn(path, selected)
        for path in ('/knowledge/findings/2', '/knowledge/findings/3', '/knowledge/findings/4', '/knowledge/genoset_audit/1'):
            self.assertNotIn(path, selected)

    def test_large_scoped_record_splits_without_losing_fields_or_expanding_subjects(self):
        r = fixture();r['knowledge']['findings'][0]['source_text'] = 'Original source passage. ' * 3000
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(r, tmp)['bundle'];m=json.loads((folder/'manifest.json').read_text())
            view=m['focus']['/conclusions/cards/0']
            self.assertGreater(len(view['files']), 1)
            closed=links(folder,view['map'])
            self.assertTrue(all('focus' in str(p) for p in closed))
            rows=[row for f in view['files'] for row in expand_part(part_data(folder/f))]
            pieces=sorted((row for row in rows if row['kind']=='json-fragment'),key=lambda x:x['offset'])
            self.assertTrue(pieces)
            self.assertEqual(json.loads(''.join(row['text'] for row in pieces)),r['knowledge']['findings'][0]['source_text'])


if __name__ == '__main__':unittest.main()
