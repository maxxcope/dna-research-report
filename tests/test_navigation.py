import hashlib
import json
from pathlib import Path
import re
import tempfile
import unittest

from genetics_app.bundle import compact
from genetics_app.navigation import Maps, MAP_BYTES, marker_links, families_for
from genetics_app.bundle import records
from genetics_app.report import export_report
from tests.test_bundle import example, rebuild


def links(folder, start):
    """Read only exported Markdown links, independently of the navigation builder."""
    seen = set()
    stack = [folder / start]
    while stack:
        path = stack.pop().resolve()
        if path in seen:
            continue
        assert path.is_relative_to(folder.resolve()), path
        assert path.is_file(), path
        seen.add(path)
        if path.suffix != '.md':
            continue
        # Read exporter-authored prose only, not links embedded in quoted source JSON.
        text = re.split(r'^(?:`{3,}|~{3,})json\n', path.read_text(), flags=re.M)[0]
        for target in re.findall(r'\]\(([^)]+)\)', text):
            if '://' not in target:
                stack.append(path.parent / target)
    return seen


class NavigationTests(unittest.TestCase):
    def test_all_maps_bounded_deep_navigation_and_source_labels_escaped(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / 'data.md').write_text('A finding')
            maps = Maps(folder)
            labels = ['Long source ' + str(i) + ' [unsafe](https://example.invalid) ' + 'x' * 150 for i in range(1800)]
            maps.tree('maps/index.md', 'Example', [(x, 'data.md') for x in labels])
            reached = links(folder, 'maps/index.md')
            self.assertIn((folder / 'data.md').resolve(), reached)
            self.assertEqual(len(reached), len(maps.inventory) + 1)
            for item in maps.inventory:
                path = folder / item['file']
                self.assertLessEqual(path.stat().st_size, MAP_BYTES)
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), item['sha256'])
                self.assertNotIn('[unsafe](https://', path.read_text())

    def test_exact_marker_links_not_substrings_or_arbitrary_prose(self):
        self.assertEqual(marker_links({'summary': 'rs123', 'id': 'rs12'}), {'rs12'})
        self.assertEqual(marker_links({'variants': [{'id': 'rs123'}], 'rule': 'and(rs9(A;A))'}), {'rs123', 'rs9'})
        report = example()
        report['knowledge'] = {'findings': [
            {'id': 'rs12', 'gene': 'SAME', 'summary': 'Example one'},
            {'id': 'rs123', 'gene': 'SAME', 'summary': 'Example two'}]}
        families, _ = families_for(report, records(report))
        self.assertNotEqual(families['rs12']['rows'], families['rs123']['rows'])
        self.assertEqual(families['rs12']['markers'], {'rs12'})

    def test_source_warning_and_clinical_record_share_marker_packet(self):
        report = example()
        report['knowledge'] = {'findings': [{'id': 'rs123', 'kind': 'variant', 'summary': 'Example finding',
                                           'assay_warnings': [{'warning': 'Keep this caution'}]}],
                               'unresolved': []}
        report['clinvar'] = [{'rsid': 'rs123', 'classification': 'conflicting', 'conditions': 'Example'}]
        report['finding_audit'] = {'flagged_findings': [{'id': 'rs123', 'flags': ['example_warning']}]}
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(report, tmp)['bundle']
            rebuilt, manifest, _ = rebuild(folder)
            self.assertEqual(compact(report), compact(rebuilt))
            matching = [p for p in manifest['parts'] if 'Example finding' in (folder / p['file']).read_text()]
            self.assertEqual(len(matching), 1)
            text = (folder / matching[0]['file']).read_text()
            self.assertIn('/finding_audit/flagged_findings/0', text)
            self.assertIn('/clinvar/0', text)
            self.assertIn('Keep this caution', text)
            self.assertNotIn('Curated trait records', (folder / 'overview.md').read_text())
            self.assertNotIn('Curated traits', (folder / 'index.md').read_text())

    def test_split_record_alias_warning_and_genoset_companions_are_reachable(self):
        report = example()
        report['knowledge'] = {'findings': [
            {'id': 'rs123', 'kind': 'variant', 'raw_ids': ['i123'], 'summary': 'A long original record',
             'source_text': 'A long source record with qualifiers. ' * 4000},
            {'id': 'gs456', 'kind': 'genoset', 'dependencies': ['rs123'], 'summary': 'Example pattern'}],
            'unresolved': [], 'genoset_audit': [{'id': 'gs456', 'rule': 'rs123(A;A)', 'reason': 'conditional'}]}
        report['finding_audit'] = {'flagged_findings': [{'id': 'i123', 'flags': ['probe warning']}],
                                  'unresolved_variants': [], 'genosets': report['knowledge']['genoset_audit']}
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(report, tmp)['bundle']
            rebuilt, manifest, _ = rebuild(folder)
            self.assertEqual(compact(report), compact(rebuilt))
            reached = links(folder, 'index.md')
            for item in manifest['parts']:
                self.assertLessEqual(item['bytes'], manifest['max_part_bytes'])
                self.assertIn((folder / item['file']).resolve(), reached)
            for item in manifest['maps']:
                self.assertLessEqual(item['bytes'], MAP_BYTES)
            seed = next(p['file'] for p in manifest['parts'] if 'Example pattern' in (folder / p['file']).read_text())
            companions = links(folder, seed)
            text = '\n'.join(p.read_text() for p in companions if p.suffix == '.md')
            self.assertIn('probe warning', text)
            self.assertIn('json-fragment', text)
            self.assertIn('/knowledge/genoset_audit/0', text)

    def test_empty_report_and_unknown_data_remain_navigable(self):
        report = example()
        report['traits'] = []
        report['unknown_field'] = {'text': 'Retain this future data'}
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(report, tmp)['bundle']
            rebuilt, manifest, _ = rebuild(folder)
            self.assertEqual(report, rebuilt)
            reached = links(folder, 'index.md')
            self.assertTrue((folder / 'maps/other/index.md').resolve() in reached)
            self.assertTrue(any('groups/other/part-' in str(p) for p in reached))


if __name__ == '__main__':
    unittest.main()
