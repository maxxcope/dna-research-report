import json
from pathlib import Path
import tempfile
import unittest

from genetics_app.bundle import at_path, expand_part
from genetics_app.report import export_report
from genetics_app.research_navigation import selections
from tests.test_bundle import example, part_data, rebuild
from tests.test_navigation import links


class ResearchNavigationTests(unittest.TestCase):
    def report(self):
        report = example()
        report['research'] = {
            'ancestry': [{'name': 'INVENTED panel', 'status': 'completed',
                          'nearest_30_population_counts': {'INVENTED': 30},
                          'points': [{'x': i, 'y': -i} for i in range(2000)],
                          'explanation': 'Similarity, not ancestry percentages.',
                          'validation': {'population_accuracy': 0.5},
                          'future_warning': 'Preserve an unknown qualification.'}],
            'maternal_lineage': {'status': 'completed', 'called_positions': 20,
                                 'results': [{'Haplogroup': 'INVENTED', 'Quality': 0.5,
                                              'Range': '1-20;' * 3000, 'Not_Found_Polys': '20G',
                                              'future_warning': 'Keep this candidate caution.'}],
                                 'explanation': 'Several branches may be compatible.'},
            'paternal_lineage': {'status': 'insufficient_data', 'candidates': [],
                                 'explanation': 'Not enough observed Y positions.'},
            'frequencies': {'status': 'completed', 'rows': [{'id': 'rs9001', 'fraction': 0.1}],
                            'explanation': 'Observed reference proportions only.'}}
        return report

    def test_summary_keeps_all_non_bulk_fields_and_original_values(self):
        report = self.report()
        views = {v['seed']: v for v in selections(report)}
        view = views['/research/ancestry/0']
        paths = view['records']
        for key in report['research']['ancestry'][0]:
            self.assertEqual('/research/ancestry/0/' + key in paths, key != 'points')
        for path, value in paths.items():
            self.assertEqual(value, at_path(report, path))
        self.assertEqual(view['details'][0][0], '/research/ancestry/0/points')

    def test_candidate_cautions_retained_while_repeated_coverage_is_optional(self):
        view = next(v for v in selections(self.report()) if v['seed'] == '/research/maternal_lineage')
        self.assertNotIn('/research/maternal_lineage/results/0/Range', view['records'])
        for field in ('Haplogroup', 'Quality', 'Not_Found_Polys', 'future_warning'):
            self.assertIn('/research/maternal_lineage/results/0/' + field, view['records'])
        self.assertEqual(view['details'][0][0], '/research/maternal_lineage/results/0/Range')

    def test_all_fields_reconstruct_and_result_plans_are_finite(self):
        report = self.report()
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(report, tmp)['bundle']
            rebuilt, manifest, _ = rebuild(folder)
            self.assertEqual(rebuilt, report)
            view = manifest['focus']['/research/ancestry/0']
            reached = links(folder, view['map'])
            self.assertEqual(reached, {(folder / f).resolve() for f in [view['map'], *view['files']]})
            data = [row for f in view['files'] for row in expand_part(part_data(folder / f))]
            self.assertFalse(any(row['path'].endswith('/points') for row in data))
            self.assertIn('Similar', '\n'.join((folder / f).read_text() for f in view['files']))
            area = (folder / 'maps/ancestry/index.md').read_text()
            self.assertIn('INVENTED panel', area)
            self.assertIn('(optional)', area)
            reached = links(folder, 'index.md')
            self.assertTrue(all((folder / x['file']).resolve() in reached for x in manifest['parts']))
            for item in manifest['maps']:
                text = (folder / item['file']).read_text()
                self.assertNotIn('- [/research/', text)

    def test_incomplete_status_cannot_be_replaced_with_success(self):
        view = next(v for v in selections(self.report()) if v['seed'] == '/research/paternal_lineage')
        self.assertEqual(view['records']['/research/paternal_lineage/status'], 'insufficient_data')
        self.assertEqual(view['records']['/research/paternal_lineage/candidates'], [])
        self.assertIn('/research/paternal_lineage/explanation', view['records'])

    def test_no_research_and_unknown_fields_are_not_invented_or_dropped(self):
        self.assertEqual(selections(example()), [])
        report = self.report()
        report['research']['maternal_lineage']['results'].append(None)
        view = next(v for v in selections(report) if v['seed'] == '/research/maternal_lineage')
        self.assertIsNone(view['records']['/research/maternal_lineage/results/1'])


if __name__ == '__main__':
    unittest.main()
