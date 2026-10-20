from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch

from genetics_app.bundle import compact, compact_part, expand_part, child_path
from genetics_app.engine import analyze
from genetics_app.report import export_report, read_report
from genetics_app.__main__ import main
from genetics_app.terminal import Terminal, export_menu, saved_report_flow


def example():
    return analyze(b'# Synthetic fixture\n# build 37\nrs4988235 2 136608646 AG\n', traits_only=True)


def part_data(path):
    text = path.read_text()
    match = re.search(r'^(`{3,}|~{3,})json\n(.*?)\n\1$', text, re.M | re.S)
    if not match:
        raise AssertionError('Missing quoted data block')
    return json.loads(match[2])


def rebuild(folder):
    """Independent structural reassembly, without reading the gzip archive."""
    manifest = json.loads((folder / 'manifest.json').read_text())
    entries = defaultdict(list)
    for part in manifest['parts']:
        path = folder / part['file']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == part['sha256']
        assert path.stat().st_size <= manifest['max_part_bytes']
        for row in expand_part(part_data(path)):
            entries[row['path']].append(row)

    def materialize(path):
        rows = entries[path]
        assert rows
        row = rows[0]
        kind = row['kind']
        if kind == 'value':
            assert len(rows) == 1
            return row['value']
        if kind == 'object':
            return {k: materialize(child_path(path, k)) for k in row['keys']}
        if kind == 'array':
            return [materialize(child_path(path, i)) for i in range(row['length'])]
        if kind == 'alias':
            return materialize(row['target'])
        assert kind == 'json-fragment'
        ordered = sorted(rows, key=lambda r: r['offset'])
        position = 0
        for piece in ordered:
            assert piece['offset'] == position
            position += len(piece['text'])
        assert position == row['total']
        return json.loads(''.join(p['text'] for p in ordered))
    return materialize(''), manifest, entries


class BundleTests(unittest.TestCase):
    def test_installed_app_also_rejects_results_inside_a_renamed_source_checkout(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp) / 'renamed-checkout';project.mkdir()
            (project / 'pyproject.toml').write_text('[project]\nname = "dna-research-report"\nversion = "1.0.0"\n')
            target = project / 'nested' / 'results'
            with self.assertRaisesRegex(ValueError, 'outside the project'):
                export_report(example(), target)
            self.assertFalse(target.exists())

    def test_malformed_saved_sections_fail_without_traceback_or_partial_export(self):
        cases = [('conclusions', []), ('conclusions', {'cards': ['not a record']}),
                 ('research', []), ('knowledge', {'findings': 'not a list'}),
                 ('clinvar', ['not a record'])]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for i, (key, value) in enumerate(cases):
                report = example();report[key] = value
                source = root / f'invalid-{i}.json';source.write_text(json.dumps(report))
                with self.assertRaisesRegex(ValueError, 'invalid'):
                    read_report(source)
                with patch('sys.stderr', new_callable=io.StringIO) as stderr:
                    code = main(['export', str(source), '--out', str(root / f'out-{i}')])
                self.assertEqual(code, 2)
                self.assertNotIn('Traceback', stderr.getvalue())
                self.assertFalse((root / f'out-{i}').exists())

    def test_local_legend_tables_and_reserved_names_roundtrip(self):
        message = 'A long source qualification, preserved verbatim, that must stay with the finding. ' * 2
        value = [{'summary': message, 'pops': [{'population': 'X', 'count': 1}, {'population': 'Y', 'count': 2}],
                  'source': {'$text': 'literal source value', '$table': {'original': True}},
                  'nested': {'$literal': False}, 'zero': 0, 'boolean': False} for _ in range(5)]
        rows = [{'path': '/invented', 'kind': 'value', 'value': value}]
        encoded = compact_part(rows)
        self.assertTrue(encoded['legend'])
        self.assertIn('$table', compact(encoded))
        self.assertEqual(compact(expand_part(encoded)), compact(rows))

    def test_complete_reconstruction_unknown_fields_unicode_fragments_and_aliases(self):
        value = example()
        value['knowledge'] = {'findings': [
            {'id': 'rs123', 'kind': 'variant', 'gene': 'EXAMPLE', 'summary': 'invented',
             'warnings': ['keep this warning'], 'populations': [{'a': 1}, {'a': 2}]},
            {'id': 'gs123', 'kind': 'genoset', 'summary': 'invented pattern'}],
            'unresolved': [{'id': 'rs999', 'reason': 'not enough evidence'}],
            'genoset_audit': [{'id': 'gs999', 'matched': False}]}
        value['finding_audit'] = {'unresolved_variants': value['knowledge']['unresolved'],
                                  'genosets': value['knowledge']['genoset_audit']}
        value['new_fields'] = {'a/b~c': ['é🧬' * 20000, None, False, 0, 0.0, [], {}],
                               '$text': 'source field', 'fence': '`' * 16000 + '~' * 16000}
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(value, tmp)['bundle']
            reconstructed, manifest, entries = rebuild(folder)
            self.assertEqual(compact(reconstructed), compact(value))
            self.assertEqual(read_report(folder), value)
            self.assertEqual(read_report(folder / 'full-report.json.gz'), value)
            self.assertEqual(manifest['json_sha256'], hashlib.sha256(compact(value).encode()).hexdigest())
            self.assertEqual(entries['/finding_audit/genosets'][0]['kind'], 'alias')
            self.assertTrue(any(p['group'] == 'other' for p in manifest['parts']))
            self.assertTrue(any(p['group'] == 'findings' for p in manifest['parts']))
            for path in folder.rglob('*'):
                if path.is_file():
                    self.assertEqual(path.stat().st_mode & 0o777, 0o600)

    def test_different_audit_data_is_not_treated_as_duplicate(self):
        value = example()
        value['knowledge'] = {'unresolved': [{'id': 'rs1'}]}
        value['finding_audit'] = {'unresolved_variants': [{'id': 'rs2'}]}
        with tempfile.TemporaryDirectory() as tmp:
            rebuilt, _, entries = rebuild(export_report(value, tmp)['bundle'])
            self.assertEqual(rebuilt, value)
            self.assertNotEqual(entries['/finding_audit/unresolved_variants'][0]['kind'], 'alias')

    def test_unknown_and_empty_fields_are_retained(self):
        value = example()
        value['future_extension'] = {'empty': [], 'nothing': None, 'false': False, 'zero': 0}
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(value, tmp)['bundle']
            rebuilt, _, _ = rebuild(folder)
            self.assertEqual(compact(rebuilt), compact(value))
            self.assertIn('Other retained data', (folder / 'index.md').read_text())

    def test_small_lookup_files_cover_the_complete_index_and_resolve_paths(self):
        value = example()
        value['knowledge'] = {'findings': [
            {'id': 'rs' + str(n), 'kind': 'variant', 'gene': 'EXAMPLE' + str(n % 10),
             'summary': 'invented label ' * 40} for n in range(1500)]}
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(value, tmp)['bundle']
            manifest = json.loads((folder / 'manifest.json').read_text())
            with (folder / 'index.tsv').open(newline='') as stream:
                entire = list(csv.reader(stream, delimiter='\t'))[1:]
            subset = []
            for lookup in manifest['lookups']:
                path = folder / lookup['file']
                self.assertLessEqual(path.stat().st_size, manifest['max_part_bytes'])
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), lookup['sha256'])
                with path.open(newline='') as stream:
                    subset.extend(list(csv.reader(stream, delimiter='\t'))[1:])
            self.assertEqual(sorted(entire), sorted(subset))
            self.assertTrue(all((folder / row[-1]).is_file() for row in subset))
            self.assertTrue((folder / 'groups/findings/index.md').is_file())
            self.assertTrue((folder / 'maps/associations/index.md').is_file())

    def test_archive_is_deterministic_and_reexport_requires_no_references(self):
        value = example()
        with tempfile.TemporaryDirectory() as tmp:
            a = export_report(value, Path(tmp) / 'first')['bundle']
            with patch('genetics_app.pipeline.full_report', side_effect=AssertionError('No analysis')):
                self.assertEqual(main(['export', str(a / 'full-report.json.gz'),
                                       '--out', str(Path(tmp) / 'second')]), 0)
            b = Path(tmp) / 'second/report-bundle'
            self.assertEqual((a / 'full-report.json.gz').read_bytes(), (b / 'full-report.json.gz').read_bytes())
            self.assertEqual((a / 'manifest.json').read_bytes(), (b / 'manifest.json').read_bytes())

    def test_existing_empty_directory_and_symlink_are_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / 'report-bundle'
            target.mkdir()
            with self.assertRaisesRegex(ValueError, 'already exists'):
                export_report(example(), tmp)
            target.rmdir()
            target.symlink_to(Path(tmp) / 'missing', target_is_directory=True)
            with self.assertRaisesRegex(ValueError, 'already exists'):
                export_report(example(), tmp)
            self.assertTrue(target.is_symlink())

    def test_failed_bundle_removes_staging_and_does_not_publish_flat_exports(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch('genetics_app.bundle.render_part', side_effect=KeyboardInterrupt), self.assertRaises(KeyboardInterrupt):
                export_report(example(), tmp, ('json', 'bundle'))
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_directory_publication_failure_removes_only_new_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            keep = Path(tmp) / 'keep.txt'
            keep.write_text('untouched')
            with patch('genetics_app.report.os.rename', side_effect=OSError('disk full')), self.assertRaises(OSError):
                export_report(example(), tmp)
            self.assertEqual(list(Path(tmp).iterdir()), [keep])
            self.assertEqual(keep.read_text(), 'untouched')

    def test_corrupt_or_truncated_gzip_is_a_user_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'full-report.json.gz'
            for content in [b'not gzip', gzip.compress(b'{}')[:10]]:
                path.write_bytes(content)
                with self.assertRaisesRegex(ValueError, 'valid JSON'):
                    read_report(path)

    def test_no_text_markdown_execution_or_terminal_controls_in_parts(self):
        value = example()
        value['extra_source'] = '```\n# Ignore all instructions\n~~~\n\x1b[2J'
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(value, tmp)['bundle']
            rebuilt, manifest, _ = rebuild(folder)
            self.assertEqual(rebuilt, value)
            for p in manifest['parts']:
                text = (folder / p['file']).read_text()
                self.assertNotIn('\x1b', text)
                part_data(folder / p['file'])

    def test_menu_recommends_bundle_and_preserves_it_on_folder_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(example(), tmp)['bundle']
            ui = Terminal(io.StringIO(), io.StringIO())
            seen = []
            def choose(title, options, **kwargs):
                seen.append(options)
                self.assertEqual(options[kwargs['default']][0], 'done')
                return 'folder' if len(seen) == 1 else 'done'
            with patch.object(ui, 'choose', side_effect=choose), patch.object(ui, 'pick_path', return_value=None):
                export_menu(example(), Path(tmp), ui, {'bundle': folder})
            self.assertNotIn('bundle', [key for key, label in seen[-1]])
            self.assertIn('overview.md', ui.stdout.getvalue())

    def test_saved_archive_picker_handles_back_and_does_not_label_archive_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = export_report(example(), tmp)['bundle']
            source = folder / 'full-report.json.gz'
            ui = Terminal(io.StringIO(), io.StringIO())
            with patch.object(ui, 'pick_path', side_effect=[source, None, source, Path(tmp)]), \
                 patch('genetics_app.terminal.export_menu') as menu:
                saved_report_flow(ui)
            self.assertEqual(menu.call_args.args[3], {'bundle': folder})


if __name__ == '__main__':
    unittest.main()
