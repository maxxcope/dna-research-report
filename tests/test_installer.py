"""Failure cases for first launch; fixtures contain no genetic observations."""
from contextlib import contextmanager
import io
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import urllib.error

from genetics_app import downloads as d, installer as i


class Response(io.BytesIO):
    def __init__(self, body, status=200, **headers):
        super().__init__(body)
        self.status = status
        self.headers = headers


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'file'
        self.quiet = lambda *args: None

    def fetch(self, responses, **kwargs):
        with patch.object(d.urllib.request, 'urlopen', side_effect=responses) as opened, patch.object(d.time, 'sleep'):
            d.download('https://example.test/data', self.path, progress=self.quiet, **kwargs)
        return opened

    def partial(self, content=b'old', validator='"v1"'):
        self.path.with_name('file.part').write_bytes(content)
        d.atomic_json(self.path.with_name('file.part.json'), dict(url='https://example.test/data', validator=validator))

    def test_download_atomic_and_verified(self):
        import hashlib
        self.fetch([Response(b'new', **{'Content-Length': '3'})], expected=hashlib.sha256(b'new').hexdigest())
        self.assertEqual(self.path.read_bytes(), b'new')
        self.assertFalse(self.path.with_name('file.part.json').exists())

    def test_changed_cached_file_is_not_accepted_or_overwritten(self):
        self.path.write_bytes(b'private fixture')
        with self.assertRaisesRegex(ValueError, 'checksum changed'):
            self.fetch([], expected='0' * 64)
        self.assertEqual(self.path.read_bytes(), b'private fixture')

    def test_invalid_download_never_becomes_final(self):
        with self.assertRaisesRegex(ValueError, 'checksum did not match'):
            self.fetch([Response(b'wrong')], expected='0' * 64)
        self.assertFalse(self.path.exists())
        self.assertFalse(self.path.with_name('file.part').exists())

    def test_guarded_range_resume(self):
        self.partial()
        opened = self.fetch([Response(b'data', 206, **{'ETag': '"v1"', 'Content-Range': 'bytes 3-6/7'})])
        request = opened.call_args.args[0]
        self.assertEqual(request.get_header('Range'), 'bytes=3-')
        self.assertEqual(request.get_header('If-range'), '"v1"')
        self.assertEqual(self.path.read_bytes(), b'olddata')

    def test_server_ignoring_range_restarts_instead_of_appending(self):
        self.partial()
        self.fetch([Response(b'replacement', **{'ETag': '"v2"'})])
        self.assertEqual(self.path.read_bytes(), b'replacement')

    def test_no_validator_restarts(self):
        self.partial(validator='')
        opened = self.fetch([Response(b'new')])
        self.assertIsNone(opened.call_args.args[0].get_header('Range'))
        self.assertEqual(self.path.read_bytes(), b'new')

    def test_incompatible_range_is_rejected(self):
        self.partial()
        with self.assertRaisesRegex(ValueError, 'incompatible partial'):
            self.fetch([Response(b'data', 206, **{'ETag': '"v1"', 'Content-Range': 'bytes 1-4/5'}) for _ in range(3)])
        self.assertFalse(self.path.exists())

    def test_bad_range_retries_from_zero(self):
        self.partial()
        opened = self.fetch([Response(b'bad', 206, **{'ETag': '"v2"', 'Content-Range': 'bytes 3-5/6'}), Response(b'new')])
        self.assertEqual(self.path.read_bytes(), b'new')
        self.assertIsNone(opened.call_args.args[0].get_header('Range'))

    def test_unsatisfied_range_restarts(self):
        self.partial()
        error = urllib.error.HTTPError('https://example.test', 416, 'range', {}, None)
        self.fetch([error, Response(b'new')])
        self.assertEqual(self.path.read_bytes(), b'new')

    def test_short_read_retries_with_matching_validator(self):
        opened = self.fetch([Response(b'old', **{'ETag': '"v1"', 'Content-Length': '7'}),
                             Response(b'data', 206, **{'ETag': '"v1"', 'Content-Range': 'bytes 3-6/7'})])
        self.assertEqual(opened.call_count, 2)
        self.assertEqual(self.path.read_bytes(), b'olddata')

    def test_http_404_is_not_retried(self):
        error = urllib.error.HTTPError('https://example.test', 404, 'missing', {}, None)
        with self.assertRaises(urllib.error.HTTPError):
            self.fetch([error])
        self.assertFalse(self.path.exists())

    def test_retry_budget_and_partial_retained(self):
        with self.assertRaises(urllib.error.URLError):
            self.fetch([urllib.error.URLError('offline')] * 3)
        self.assertFalse(self.path.exists())

    def test_pin_size_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError, 'size differs'):
            self.fetch([Response(b'abc', **{'Content-Length': '3'})], size=5)

    def test_disk_error_does_not_loop(self):
        with patch.object(d, 'atomic_json', side_effect=OSError(28, 'No space left')):
            with self.assertRaises(OSError):
                self.fetch([Response(b'data')])
        self.assertFalse(self.path.exists())

    def test_interrupt_keeps_final_unpublished(self):
        with self.assertRaises(KeyboardInterrupt):
            self.fetch([KeyboardInterrupt()])
        self.assertFalse(self.path.exists())

    def test_malformed_checkpoint_restarts(self):
        self.partial()
        self.path.with_name('file.part.json').write_text('{bad')
        self.fetch([Response(b'new')])
        self.assertEqual(self.path.read_bytes(), b'new')

    def test_valid_json_with_wrong_checkpoint_shape_restarts(self):
        self.partial()
        self.path.with_name('file.part.json').write_text('["broken checkpoint"]')
        self.fetch([Response(b'new')])
        self.assertEqual(self.path.read_bytes(), b'new')


class InstallerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_cache_inside_project_rejected_even_through_symlink(self):
        (self.root / 'link').symlink_to(i.ROOT, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'outside this project'):
            i.external_folder(self.root / 'link/cache')
        self.assertEqual(i.external_folder(self.root / 'outside'), (self.root / 'outside').resolve())

    def test_lock_released_after_failure_and_duplicate_rejected(self):
        with self.assertRaisesRegex(RuntimeError, 'interrupted'):
            with i.setup_lock(self.root):
                with self.assertRaisesRegex(ValueError, 'already running'):
                    with i.setup_lock(self.root):
                        self.fail('duplicate setup admitted')
                raise RuntimeError('interrupted')
        with i.setup_lock(self.root):
            pass

    def test_unsupported_platform_fails_before_install(self):
        with patch.object(i.platform, 'system', return_value='Windows'):
            with self.assertRaisesRegex(ValueError, 'supports macOS'):
                i.platform_key()

    def test_plink2_and_old_java_are_rejected(self):
        for name, text in [('plink', 'PLINK v2.0'), ('java', 'java version "1.7.0"')]:
            with patch.object(i.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, text, '')):
                self.assertFalse(i.executable_works('tool', name))

    def test_java_and_plink_version_detection(self):
        for name, text in [('plink', 'PLINK v1.9.0-b.7.13 64-bit'), ('java', 'openjdk version "21.0.12.1"')]:
            with patch.object(i.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, text, '')):
                self.assertTrue(i.executable_works('tool', name))

    def profile(self):
        data = self.root / 'public-reference'
        data.write_text('invented reference fixture')
        f = str(data)
        prefix = str(self.root / 'panel')
        files = [data]
        for suffix in ['bed', 'bim', 'fam']:
            p = Path(prefix + '.' + suffix)
            p.write_text('invented public panel fixture')
            files.append(p)
        records = [dict(file=str(p), sha256=d.sha256(p)) for p in files]
        config = dict(schema=1, knowledge=f, clinvar=f, gwas=f,
                      quality=dict(clusters=f, badalleles=f),
                      frequencies=dict(prefix=prefix, labels=f),
                      paternal=dict(tree=f, chain=f),
                      ancestry=[dict(prefix=prefix, labels=f, plink=f)],
                      maternal=dict(directory=str(self.root), java=f, tree='fixture@1'),
                      scores=[dict(id='fixture', path=f, build='37')],
                      imputation=dict(windows=[dict(file=f, map=f)], java=f, beagle=f, score=f),
                      phase=dict(window=dict(file=f, map=f), java=f, beagle=f, sites=[{}]),
                      reference_files=records)
        d.atomic_json(self.root / 'pipeline.json', config)
        return data

    @contextmanager
    def profile_mocks(self):
        with patch('genetics_app.references.metadata', return_value={'build': '37'}), patch.object(i, 'executable_works', return_value=True):
            yield

    def test_ready_profile_reuses_metadata_without_rehashing_data(self):
        self.profile()
        marker = self.root / 'ready.json'
        with self.profile_mocks():
            self.assertTrue(i.validate_profile(self.root, marker))
            original = i.sha256
            def guard(path):
                self.assertEqual(Path(path).name, 'pipeline.json')
                return original(path)
            with patch.object(i, 'sha256', side_effect=guard):
                self.assertTrue(i.validate_profile(self.root, marker))

    def test_changed_reference_never_gets_new_trusted_hash(self):
        data = self.profile()
        marker = self.root / 'ready.json'
        with self.profile_mocks():
            i.validate_profile(self.root, marker)
            before = marker.read_bytes()
            data.write_text('changed public reference')
            with self.assertRaisesRegex(ValueError, 'checksum changed'):
                i.validate_profile(self.root, marker)
            self.assertEqual(marker.read_bytes(), before)

    def test_missing_reference_requests_resume(self):
        self.profile().unlink()
        with self.profile_mocks():
            self.assertFalse(i.validate_profile(self.root, self.root / 'ready.json'))

    def test_incomplete_profile_is_not_success(self):
        d.atomic_json(self.root / 'pipeline.json', dict(schema=1))
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            i.validate_profile(self.root, self.root / 'ready.json')

    def test_malformed_profile_is_actionable(self):
        d.atomic_json(self.root / 'pipeline.json', ['broken profile fixture'])
        with self.assertRaisesRegex(ValueError, 'Malformed reference profile'):
            i.validate_profile(self.root, self.root / 'ready.json')

    def test_untracked_configured_file_is_rejected(self):
        self.profile()
        profile = self.root / 'pipeline.json'
        config = json.loads(profile.read_text())
        config['gwas'] = str(self.root / 'untracked')
        d.atomic_json(profile, config)
        with self.assertRaisesRegex(ValueError, 'missing from the hash inventory'):
            i.validate_profile(self.root, self.root / 'ready.json')

    def test_malformed_file_inventory_is_rejected(self):
        self.profile()
        profile = self.root / 'pipeline.json'
        config = json.loads(profile.read_text())
        config['reference_files'][0]['sha256'] = 'not a hash'
        d.atomic_json(profile, config)
        with self.assertRaisesRegex(ValueError, 'Malformed reference file inventory'):
            i.validate_profile(self.root, self.root / 'ready.json')

    def test_no_ready_marker_if_hash_check_fails(self):
        self.profile().write_text('bad')
        marker = self.root / 'ready.json'
        with self.profile_mocks(), self.assertRaisesRegex(ValueError, 'checksum changed'):
            i.validate_profile(self.root, marker)
        self.assertFalse(marker.exists())

    def test_tool_archive_traversal_rejected(self):
        archive = self.root / 'bad.tar.gz'
        with tarfile.open(archive, 'w:gz') as target:
            member = tarfile.TarInfo('../outside')
            member.size = 3
            target.addfile(member, io.BytesIO(b'bad'))
        dest = self.root / 'extract'
        dest.mkdir()
        with self.assertRaises((tarfile.FilterError, ValueError)):
            i.extract_tool(archive, dest, 'java')
        self.assertFalse((self.root / 'outside').exists())

    def test_tool_archive_escaping_symlink_rejected(self):
        archive = self.root / 'bad.tar.gz'
        with tarfile.open(archive, 'w:gz') as target:
            member = tarfile.TarInfo('escape')
            member.type = tarfile.SYMTYPE
            member.linkname = '../../outside'
            target.addfile(member)
        dest = self.root / 'extract'
        dest.mkdir()
        with self.assertRaises((tarfile.FilterError, ValueError)):
            i.extract_tool(archive, dest, 'java')

    def test_too_little_disk_stops_before_tool_download(self):
        with patch.object(i, 'ensure_environment', return_value=Path('python')), patch.object(i, 'validate_profile', return_value=False), patch.object(i.shutil, 'disk_usage', return_value=type('Disk', (), {'free': 1})()), patch.object(i, 'ensure_tool') as tools:
            with self.assertRaisesRegex(ValueError, '10 GiB'):
                i.setup(self.root, self.root / 'refs', 'uv', self.root / 'log')
            tools.assert_not_called()

    def test_completed_setup_does_not_download_or_prepare(self):
        with patch.object(i, 'ensure_environment', return_value=Path('python')), patch.object(i, 'validate_profile', return_value=True), patch.object(i, 'ensure_tool') as tools, patch.object(i, 'run') as run:
            self.assertEqual(i.setup(self.root, self.root / 'refs', 'uv', self.root / 'log'), Path('python'))
            tools.assert_not_called()
            run.assert_not_called()

    def test_ready_environment_skips_package_manager(self):
        with patch.object(i, 'environment_ready', return_value=True), patch.object(i, 'run') as run:
            self.assertEqual(i.ensure_environment(self.root, 'uv', self.root / 'log'), self.root / 'environment/bin/python')
            run.assert_not_called()

    def test_partial_environment_is_recreated_and_verified(self):
        (self.root / 'environment').mkdir()
        (self.root / 'environment/incomplete').write_text('interrupted install fixture')
        with patch.object(i, 'environment_ready', side_effect=[False, True]), patch.object(i, 'run') as run:
            i.ensure_environment(self.root, 'uv', self.root / 'log')
            self.assertEqual(run.call_count, 2)
            self.assertIn('--require-hashes', run.call_args.args[0])
            self.assertFalse((self.root / 'environment/incomplete').exists())

    def test_bad_tool_override_stops_before_download(self):
        with patch.dict(i.os.environ, {'DNA_RESEARCH_PLINK': '/nonexistent-fixture-tool'}), patch.object(i, 'executable_works', return_value=False), patch.object(i, 'download') as download:
            with self.assertRaisesRegex(ValueError, 'configured plink'):
                i.ensure_tool(self.root, 'plink', {})
            download.assert_not_called()

    def test_changed_tool_payload_is_rejected(self):
        spec = dict(sha256='1' * 64)
        folder = self.root / 'tools' / ('plink-' + spec['sha256'][:16])
        folder.mkdir(parents=True)
        (folder / 'plink').write_bytes(b'changed executable fixture')
        d.atomic_json(folder / 'installed.json', dict(executable='plink', archive_sha256=spec['sha256'], files=[dict(file='plink', sha256='0' * 64)]))
        with patch.object(i, 'download') as download:
            with self.assertRaisesRegex(ValueError, 'installed plink files changed'):
                i.ensure_tool(self.root, 'plink', spec)
            download.assert_not_called()

    def test_failed_setup_returns_error_without_ready_marker(self):
        with patch.dict(i.os.environ, {'DNA_RESEARCH_HOME': str(self.root)}), patch.object(i, 'platform_key', return_value='darwin-arm64'), patch.object(i, 'setup', side_effect=ValueError('offline fixture')), patch('sys.stdout', new=io.StringIO()):
            self.assertEqual(i.main(['--setup-only']), 1)
        self.assertFalse((self.root / 'references-ready.json').exists())
        with i.setup_lock(self.root):
            pass

    def test_cancelled_setup_returns_130_and_unlocks(self):
        with patch.dict(i.os.environ, {'DNA_RESEARCH_HOME': str(self.root)}), patch.object(i, 'platform_key', return_value='darwin-arm64'), patch.object(i, 'setup', side_effect=KeyboardInterrupt()), patch('sys.stderr', new=io.StringIO()):
            self.assertEqual(i.main(['--setup-only']), 130)
        with i.setup_lock(self.root):
            pass

    def test_app_only_explicitly_bypasses_reference_setup(self):
        with patch.dict(i.os.environ, {'DNA_RESEARCH_HOME': str(self.root)}), patch.object(i, 'platform_key', return_value='darwin-arm64'), patch.object(i, 'setup') as setup, patch.object(i, 'ensure_environment', return_value=Path('python')), patch.object(i, 'launch') as launch:
            i.main(['--app-only'])
            setup.assert_not_called()
            self.assertEqual(launch.call_args.args[1], [])

    def test_reference_download_rejects_mixed_source_snapshot(self):
        from genetics_app.setup import download
        path = self.root / 'public-file'
        with self.assertRaisesRegex(ValueError, 'source changed'):
            download('https://example.test/new', path, {'public-file': dict(url='https://example.test/old', sha256='0' * 64)})

    def test_shared_reference_lock_prevents_different_installers(self):
        from genetics_app.setup import prepare
        def nested(*args):
            with self.assertRaisesRegex(ValueError, 'already running for this folder'):
                prepare(self.root, 'plink', 'java')
        with patch('genetics_app.setup._prepare', side_effect=nested):
            prepare(self.root, 'plink', 'java')
        with patch('genetics_app.setup._prepare', return_value='ready'):
            self.assertEqual(prepare(self.root, 'plink', 'java'), 'ready')


class ReferenceWindowSetupTests(unittest.TestCase):
    def setUp(self):
        import pysam
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/'public-fixture.vcf.gz'
        header=pysam.VariantHeader()
        header.contigs.add('1',length=1000)
        header.formats.add('GT',1,'String','Genotype')
        header.add_sample('INVENTED_REF_A');header.add_sample('INVENTED_REF_B')
        with pysam.VariantFile(str(self.source),'wz',header=header) as output:
            rec=output.new_record(contig='1',start=99,alleles=('A','G'))
            rec.samples['INVENTED_REF_A']['GT']=(0,0)
            rec.samples['INVENTED_REF_B']['GT']=(0,1)
            output.write(rec)
        pysam.tabix_index(str(self.source),preset='vcf')
        self.window=dict(name='invented',chrom='1',start=90,end=110)
        self.samples=['INVENTED_REF_B']

    def prepare(self,manifest,checkpoint=lambda:None,samples=None):
        from genetics_app.setup import prepare_windows
        with patch('genetics_app.setup.VCF',str(self.source)), patch('genetics_app.setup.download',return_value=Path(str(self.source)+'.tbi')):
            return prepare_windows([self.window],self.root,self.samples if samples is None else samples,manifest,checkpoint)

    def test_window_has_exact_sample_subset_index_and_checkpoint(self):
        import pysam
        manifest={};calls=[]
        out=self.prepare(manifest,lambda:calls.append(dict(manifest)))
        self.assertEqual(len(calls),2)
        self.assertIn('invented',calls[-1])
        path=Path(out[0]['file']);self.assertTrue(Path(str(path)+'.tbi').is_file())
        with pysam.VariantFile(path) as source:
            self.assertEqual(list(source.header.samples),self.samples)
            records=list(source.fetch('1',89,110))
            self.assertEqual(len(records),1)
            self.assertEqual(records[0].samples[self.samples[0]]['GT'],(0,1))

    def test_checkpoint_interrupt_resumes_without_remote_source(self):
        manifest={}
        def interrupt():
            if 'invented' in manifest:raise KeyboardInterrupt()
        with self.assertRaises(KeyboardInterrupt):self.prepare(manifest,interrupt)
        self.assertIn('invented',manifest)
        self.source.unlink()
        self.assertEqual(len(self.prepare(manifest)),1)

    def test_modified_window_is_rejected_before_indexing(self):
        manifest={};out=self.prepare(manifest)
        Path(out[0]['file']).write_bytes(b'changed cached fixture')
        with self.assertRaisesRegex(ValueError,'Cached window checksum changed'):
            self.prepare(manifest)

    def test_different_sample_order_is_not_reused(self):
        manifest={};self.prepare(manifest)
        with self.assertRaisesRegex(ValueError,'unexpected sample order'):
            self.prepare(manifest,samples=['INVENTED_REF_A'])


if __name__ == '__main__':
    unittest.main()
