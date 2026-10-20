import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from contextlib import contextmanager
from unittest.mock import patch
from genetics_app.engine import analyze
from genetics_app.parser import detect_build
from genetics_app.report import export_report, report_lines, safe_text, output_folder, read_report
from genetics_app.progress import Progress
from genetics_app.terminal import Terminal, guided, path_value
from genetics_app.__main__ import main
from genetics_app import extras

DATA = b'# Synthetic fixture\n# build 37\n# Genotype calls use the forward strand.\nrs4988235 2 136608646 AG\n'


def report():
    return analyze(DATA, traits_only=True)


class ExportTests(unittest.TestCase):
    def test_project_local_results_are_rejected_before_creating_files(self):
        project=Path(__file__).resolve().parents[1]
        target=project/'private-result-regression'
        with self.assertRaisesRegex(ValueError,'outside the project'):
            export_report(report(),target)
        self.assertFalse(target.exists())

    def test_saved_report_rejects_malformed_types_without_traceback(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.json'
            for bad in [[], {'input_sha256':'x','traits':{},'clinvar':[],'gwas':[]},
                        dict(report(),pipeline={'status':'completed','stages':[None]}),
                        dict(report(),knowledge={'counts':{'variant_summaries':'bad'}})]:
                path.write_text(json.dumps(bad))
                with self.subTest(value=bad),self.assertRaises(ValueError):read_report(path)

    def test_exports_only_offer_unsaved_formats_and_finish_after_success(self):
        from genetics_app.terminal import export_menu
        with tempfile.TemporaryDirectory() as tmp:
            saved=export_report(report(),tmp,('json',))
            ui=Terminal(io.StringIO(),io.StringIO());calls=[]
            def choose(prompt, options, **kwargs):
                keys=[key for key,label in options];calls.append(keys)
                self.assertNotIn('json',keys)
                if len(calls)==1:return 'all'
                self.assertEqual(keys,['folder','done'])
                self.assertEqual(keys[kwargs['default']],'done')
                return 'done'
            with patch.object(ui,'choose',side_effect=choose):export_menu(report(),Path(tmp),ui,saved)
            self.assertTrue(Path(tmp,'report.txt').is_file())

    def test_quit_after_changing_export_folder_names_already_saved_reports(self):
        from genetics_app.terminal import export_menu,QuitFlow
        with tempfile.TemporaryDirectory() as tmp:
            saved=export_report(report(),tmp,('json',));ui=Terminal(io.StringIO(),io.StringIO())
            with patch.object(ui,'choose',side_effect=['folder',QuitFlow]),patch.object(ui,'pick_path',return_value=Path(tmp)/'another'):
                export_menu(report(),Path(tmp),ui,saved)
            self.assertIn(str(saved['json']),ui.stdout.getvalue())
            self.assertNotIn('No report was saved',ui.stdout.getvalue())

    def test_returning_to_export_folder_remembers_completed_formats(self):
        from genetics_app.terminal import export_menu
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);saved=export_report(report(),root,('json',))
            ui=Terminal(io.StringIO(),io.StringIO());calls=[]
            def choose(prompt,options,**kwargs):
                keys=[k for k,label in options];calls.append(keys)
                if len(calls)==5:
                    self.assertEqual(keys,['folder','done'])
                    return 'done'
                return ['all','folder','folder','folder'][len(calls)-1]
            with patch.object(ui,'choose',side_effect=choose),patch.object(ui,'pick_path',side_effect=[root/'other',root,root]):
                export_menu(report(),root,ui,saved)
    def test_all_fields_and_late_records_preserved_in_every_format(self):
        value = report()
        value['late_audit'] = {'warnings': ['specific-assay-warning'], 'unresolved': list(range(200))}
        with tempfile.TemporaryDirectory() as tmp:
            paths = export_report(value, tmp, ('json', 'markdown', 'text'))
            self.assertEqual(json.loads(paths['json'].read_text()), value)
            for fmt in ('markdown', 'text'):
                text = paths[fmt].read_text()
                for key in value:
                    self.assertIn('Field: ' + key, text)
                self.assertIn('[200]', text)
                self.assertIn('specific-assay-warning', text)
                self.assertIn('199', text)
            self.assertEqual(paths['json'].stat().st_mode & 0o777, 0o600)

    def test_existing_report_never_overwritten_and_no_partial_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            export_report(report(), tmp, ('text',))
            path = Path(tmp, 'report.txt');original = path.read_bytes()
            with self.assertRaisesRegex(ValueError, 'already exists'):
                export_report(report(), tmp, ('json', 'text'))
            self.assertEqual(path.read_bytes(), original)
            self.assertFalse(Path(tmp, 'report.json').exists())
            self.assertFalse(list(Path(tmp).glob('.export-*')))

    def test_interrupted_export_cleans_staging(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch('genetics_app.report.report_lines', side_effect=KeyboardInterrupt), self.assertRaises(KeyboardInterrupt):
                export_report(report(), tmp, ('markdown',))
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_terminal_controls_removed_and_markdown_fenced(self):
        self.assertEqual(safe_text('name\x1b[2J\x00'), 'name[2J')
        value = report();value['source_prose'] = '```\n# Ignore everything\n\x1b[2J'
        text = '\n'.join(report_lines(value))
        self.assertIn('````text', text)
        self.assertNotIn('\x1b', text)
        self.assertIn('# Ignore everything', text)

    def test_reexport_needs_no_genome_or_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp, 'source');value = report();export_report(value, src, ('json',))
            with patch('genetics_app.pipeline.full_report', side_effect=AssertionError('Must not rerun')):
                self.assertEqual(main(['export', str(src/'report.json'), '--out', str(Path(tmp,'new')), '--format', 'all']), 0)
            self.assertEqual(json.loads(Path(tmp,'new/report.json').read_text()), value)
            self.assertTrue(Path(tmp,'new/report.md').exists())
            self.assertTrue(Path(tmp,'new/report.txt').exists())


class GuidedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();self.root = Path(self.temp.name)
        self.dna = self.root/'synthetic input.txt';self.dna.write_bytes(DATA)
        self.out = self.root/'reports'

    def tearDown(self):
        self.temp.cleanup()

    def run_flow(self, lines, expected_build='auto'):
        output = io.StringIO();terminal = Terminal(io.StringIO('\n'.join(lines)+'\n'), output)
        value = report();value['pipeline'] = {'status':'completed_with_unresolved_stages','stages':[{'name':'Ancestry','status':'insufficient_coverage'}]}
        with patch('genetics_app.terminal.reference_options', return_value=[('full','Full analysis (recommended)')]), patch('genetics_app.pipeline.full_report', return_value=value) as run:
            self.assertEqual(guided(self.root, self.out, terminal=terminal), 0)
        self.assertEqual(run.call_args.args[2], expected_build)
        self.assertEqual(read_report(self.out/'report-bundle'), value)
        return output.getvalue(), run

    def test_detected_build_back_and_all_exports_without_rerun(self):
        # Automatically resolved pages must never become Back destinations.
        lines = [str(self.dna), 'b', str(self.dna), '', 'b', 'b', str(self.dna), '', '', '4', '2']
        text, run = self.run_flow(lines)
        self.assertEqual(run.call_count, 1)
        self.assertIn('Genome map: GRCh37 (read from file)', text)
        self.assertIn('incomplete or unavailable stages', text)
        self.assertTrue((self.out/'report.md').exists());self.assertTrue((self.out/'report.txt').exists())

    def test_declared_build_and_full_references_require_no_confirmation(self):
        text, _ = self.run_flow([str(self.dna), '', '', '5'])
        self.assertIn('Genome map: GRCh37 (read from file)', text)
        self.assertNotIn('Which map does your provider specify?', text)
        self.assertNotIn('Available analysis', text)

    def test_detected_format_skips_confirmation_and_goes_to_output(self):
        text, run = self.run_flow([str(self.dna), '', '', '5'])
        self.assertEqual(run.call_count, 1)
        self.assertNotIn('DNA letter format', text)
        self.assertNotIn('Continue with this file', text)
        self.assertNotIn('Check your file format', text)
        self.assertIn('Save reports in', text)

    def test_reverse_header_returns_to_file_before_analysis(self):
        bad=self.root/'reverse.txt';bad.write_bytes(DATA.replace(b'forward strand',b'reverse strand'))
        text, run=self.run_flow([str(bad), str(self.dna), '', '', '5'])
        self.assertEqual(run.call_count,1)
        self.assertIn('Get a forward/plus export',text)

    def test_export_failure_keeps_result_and_retries_without_analysis(self):
        real_export = export_report
        failed = []
        def save(value, folder, formats):
            if formats == ('markdown',) and not failed:
                failed.append(True)
                raise OSError('Simulated full disk')
            return real_export(value, folder, formats)
        with patch('genetics_app.terminal.export_report', side_effect=save):
            text, run = self.run_flow([str(self.dna), '', '', '2', '4', '2'])
        self.assertEqual(run.call_count, 1)
        self.assertIn('completed analysis is still available', text)
        self.assertTrue((self.out/'report.md').exists())

    def test_missing_build_requires_explicit_provider_choice(self):
        self.dna.write_bytes(DATA.replace(b'# build 37\n', b''))
        text, _ = self.run_flow([str(self.dna), '1', '', '', '5'], '37')
        self.assertIn('Your file does not name its genome build', text)

    def test_conflicting_header_returns_to_file(self):
        bad = self.root/'bad.txt';bad.write_bytes(DATA + b'# GRCh38\n')
        text, _ = self.run_flow([str(bad), str(self.dna), '', '', '5'])
        self.assertIn('Conflicting genome builds', text)

    def test_early_quit_does_not_create_report(self):
        terminal = Terminal(io.StringIO('q\n'), io.StringIO())
        from genetics_app.terminal import QuitFlow
        with self.assertRaises(QuitFlow):guided(self.root, self.out, terminal=terminal)
        self.assertFalse(self.out.exists())

    def test_missing_references_cannot_silently_run_traits(self):
        terminal = Terminal(io.StringIO(str(self.dna)+'\nq\n'), io.StringIO())
        from genetics_app.terminal import QuitFlow
        with self.assertRaises(QuitFlow), patch('genetics_app.pipeline.full_report') as run:
            guided(self.root, self.out, terminal=terminal)
        run.assert_not_called();self.assertFalse(self.out.exists())

    def test_review_edits_can_be_cancelled_without_losing_selections(self):
        ui=Terminal(io.StringIO(),io.StringIO())
        with patch.object(ui,'pick_path',side_effect=[self.dna,self.out,None,None,None]),patch.object(ui,'choose',side_effect=['output','file','refs','run','done']),patch('genetics_app.terminal.reference_options',return_value=[('full','Full')]),patch('genetics_app.pipeline.full_report',return_value=report()) as run:
            self.assertEqual(guided(self.root,self.out,terminal=ui),0)
        self.assertEqual(run.call_count,1);self.assertEqual(run.call_args.args[0],DATA)
        self.assertTrue((self.out/'report-bundle/full-report.json.gz').exists())

    def test_quit_before_generate_does_not_create_report_directory(self):
        from genetics_app.terminal import QuitFlow
        ui=Terminal(io.StringIO(str(self.dna)+'\n\nq\n'),io.StringIO())
        with patch('genetics_app.terminal.reference_options',return_value=[('full','Full')]),self.assertRaises(QuitFlow):
            guided(self.root,self.out,terminal=ui)
        self.assertFalse(self.out.exists())

    def test_analysis_failure_can_retry_from_review(self):
        ui=Terminal(io.StringIO(str(self.dna)+'\n\n\n\n5\n'),io.StringIO())
        with patch('genetics_app.terminal.reference_options',return_value=[('full','Full')]),patch('genetics_app.pipeline.full_report',side_effect=[ValueError('Checksum changed'),report()]) as run:
            self.assertEqual(guided(self.root,self.out,terminal=ui),0)
        self.assertEqual(run.call_count,2)
        self.assertIn('Analysis failed: Checksum changed',ui.stdout.getvalue())

    def test_unknown_build_back_returns_to_the_build_question(self):
        self.dna.write_bytes(DATA.replace(b'# build 37\n',b''))
        text,_=self.run_flow([str(self.dna),'1','b','2','','','5'],'38')
        self.assertEqual(text.count('Which map does your provider specify?'),2)

    def test_reference_fallback_never_silently_selects_limited_analysis(self):
        from genetics_app.terminal import select_references
        ui=Terminal(io.StringIO(),io.StringIO())
        with patch('genetics_app.terminal.reference_options',return_value=[('basic','Limited')]),patch.object(ui,'choose',return_value='back') as choose:
            self.assertIsNone(select_references(self.root,'38',ui))
        choose.assert_called_once()

    def test_confirmed_full_profile_is_automatic_even_when_basic_is_available(self):
        from genetics_app.terminal import select_references
        ui=Terminal(io.StringIO(),io.StringIO())
        with patch('genetics_app.terminal.reference_options',return_value=[('full','Full'),('basic','Limited')]),patch.object(ui,'choose',side_effect=AssertionError('No decision needed')):
            self.assertEqual(select_references(self.root,'37',ui),(self.root,'full',False))

    def test_path_with_spaces_shell_quotes_and_drag_escape(self):
        self.assertEqual(path_value(str(self.dna)), self.dna.resolve())
        self.assertEqual(path_value("'"+str(self.dna)+"'"), self.dna.resolve())
        self.assertEqual(path_value(str(self.dna).replace(' ', '\\ ')), self.dna.resolve())

    def test_missing_and_conflicting_detection(self):
        self.assertIsNone(detect_build(b'rs1 1 100 AA'))
        self.assertEqual(detect_build(b'# hg19\nrs1 1 100 AA'), '37')
        with self.assertRaises(ValueError):detect_build(b'# build 36')
        with self.assertRaises(ValueError):detect_build(b'# hg19 and hg38')


class ProgressAndCancellationTests(unittest.TestCase):
    def test_machine_events_are_actual_counts_without_ansi(self):
        out = io.StringIO()
        with Progress(out, mode='json') as p:
            p.event(dict(kind='stage_started', name='Stage', completed=0, total=2))
            p.event(dict(kind='detail', message='Imputation windows 3/7'))
            p.event(dict(kind='stage_finished', name='Stage', completed=1, total=2, status='completed'))
        events = [json.loads(s) for s in out.getvalue().splitlines()]
        self.assertEqual(events[0]['completed'], 0);self.assertEqual(events[-1]['completed'], 1)
        self.assertIn('3/7', events[1]['message']);self.assertNotIn('\x1b', out.getvalue())

    def test_tool_cancellation_reaps_process_and_allows_next_run(self):
        extras.reset_commands();failures = []
        def worker():
            try:extras.run_command([sys.executable, '-c', 'import time;time.sleep(30)'], '.', timeout=60)
            except ValueError as e:failures.append(str(e))
        t = threading.Thread(target=worker);t.start()
        for _ in range(100):
            with extras._command_lock:active = list(extras._commands)
            if active:break
            time.sleep(.01)
        self.assertTrue(active);extras.cancel_commands();t.join(5)
        self.assertFalse(t.is_alive());self.assertIsNotNone(active[0].poll());self.assertTrue(failures)
        extras.reset_commands()
        self.assertEqual(extras.run_command([sys.executable, '-c', 'print("ready")'], '.').strip(), 'ready')

    def test_no_arguments_noninteractive_prints_help(self):
        with patch('sys.stdin', io.StringIO()), patch('sys.stdout', io.StringIO()):
            self.assertEqual(main([]), 2)


class TerminalRestorationTests(unittest.TestCase):
    def test_new_screen_scrolls_previous_page_into_history(self):
        output = io.StringIO();ui = Terminal(io.StringIO(), output);ui.arrows = True
        ui.say('Previous step')
        with patch('genetics_app.terminal.shutil.get_terminal_size', return_value=os.terminal_size((80, 24))):
            ui.screen()
        ui.say('Next step')
        self.assertEqual(output.getvalue(), 'Previous step\n' + '\r\n'*24 + '\x1b[HNext step\n')
        for erase in ('\x1b[2J', '\x1b[3J', '\x1b[?1049h', '\x1bc'):
            self.assertNotIn(erase, output.getvalue())

    def test_redirected_output_has_no_screen_controls(self):
        output = io.StringIO();ui = Terminal(io.StringIO(), output)
        ui.say('Previous step');ui.screen();ui.say('Next step')
        self.assertEqual(output.getvalue(), 'Previous step\nNext step\n')

    def test_export_error_stays_visible_after_new_screen(self):
        from genetics_app.terminal import export_menu
        output=io.StringIO();ui=Terminal(io.StringIO(),output)
        def clear():
            output.seek(0);output.truncate()
        def choose(*args, **kwargs):
            self.assertIn('Could not save JSON', output.getvalue())
            self.assertIn('Choose another export folder', output.getvalue())
            return 'done'
        with patch.object(ui, 'screen', side_effect=clear),patch.object(ui, 'choose', side_effect=choose):
            export_menu(report(), Path('/tmp/invented-output'), ui, {}, notice='Could not save JSON\nChoose another export folder')

    def test_progress_cursor_restored_on_interrupt(self):
        class TTY(io.StringIO):
            def isatty(self):return True
        output = TTY()
        with self.assertRaises(KeyboardInterrupt):
            with Progress(output) as p:
                p.event(dict(kind='stage_started', name='Actual stage', completed=2, total=12))
                raise KeyboardInterrupt()
        self.assertIn('\x1b[?25l', output.getvalue())
        self.assertIn('\x1b[?25h', output.getvalue())
        self.assertNotIn('12/12 stages processed', output.getvalue())


@unittest.skipIf(os.name == 'nt', 'POSIX terminal integration')
class PasteInputTests(unittest.TestCase):
    @contextmanager
    def session(self, expression="ui.line('Path')", setup=''):
        import pty
        import termios
        master, slave = pty.openpty()
        previous = termios.tcgetattr(slave)
        script = ('import json,signal\nsignal.signal(signal.SIGINT,signal.default_int_handler)\n'
                  'signal.pthread_sigmask(signal.SIG_UNBLOCK,{signal.SIGINT})\n'
                  'from genetics_app.terminal import Terminal\n'
                  + setup + '\n'
                  'ui=Terminal()\nvalue=' + expression + '\n'
                  'print("RESULT="+json.dumps(value),flush=True)\n')
        process = subprocess.Popen([sys.executable, '-u', '-c', script], stdin=slave,
                                   stdout=slave, stderr=slave, cwd=Path(__file__).resolve().parents[1])
        self.received = b''
        self.master = master
        try:
            yield process, slave, previous
        finally:
            if process.poll() is None:
                process.kill()
            process.wait(timeout=5)
            os.close(master)
            os.close(slave)

    def receive(self, needle=None, timeout=3):
        import select
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if needle is not None and needle in self.received:
                return self.received
            if not select.select([self.master], [], [], min(.05, max(0, deadline-time.monotonic())))[0]:
                continue
            try:
                chunk = os.read(self.master, 65536)
            except OSError:
                break
            if not chunk:
                break
            self.received += chunk
        if needle is not None:
            self.assertIn(needle, self.received)
        return self.received

    def send(self, text):
        os.write(self.master, text.encode('utf-8'))

    def result(self, expected, process):
        self.receive(('RESULT='+json.dumps(expected)).encode())
        self.assertEqual(process.wait(timeout=3), 0)

    def assert_restored(self, slave, previous):
        import termios
        current = termios.tcgetattr(slave)
        # macOS may set its pending-input bookkeeping flag after restoration.
        current[3] &= ~getattr(termios, 'PENDIN', 0)
        previous[3] &= ~getattr(termios, 'PENDIN', 0)
        self.assertEqual(current, previous)

    def test_pasted_trailing_newline_waits_for_enter_and_can_be_extended(self):
        with self.session() as (process, slave, previous):
            self.receive(b'Path: ')
            self.send('\x1b[200~/tmp/invented\n\x1b[201~')
            self.receive(b'/tmp/invented')
            self.receive(timeout=.15)
            self.assertIsNone(process.poll())
            self.assertNotIn(b'RESULT=', self.received)
            self.send(' sample.txt\r')
            self.result('/tmp/invented sample.txt', process)
            self.assert_restored(slave, previous)
            self.assertIn(b'\x1b[?2004l', self.received)

    def test_fragmented_multiline_paste_does_not_submit_or_run_shortcuts(self):
        with self.session() as (process, _, __):
            self.receive(b'Path: ')
            self.send('\x1b[200~/tmp/')
            time.sleep(.03)
            self.send('q\r\nb\r\n')
            time.sleep(.03)
            self.send('sample.txt\r\n\x1b[201~')
            self.receive(timeout=.15)
            self.assertIsNone(process.poll())
            self.send('\r')
            self.result('/tmp/q b sample.txt', process)

    def test_unicode_dragged_path_and_cursor_editing(self):
        with self.session() as (process, _, __):
            self.receive(b'Path: ')
            self.send('\x1b[200~/tmp/café\\ sample.tx\n\x1b[201~')
            self.send('X\x7ft\x1b[D\x1b[3~T\r')
            self.result('/tmp/café\\ sample.txT', process)

    def test_paste_over_menu_cannot_select_quit_back_or_enter(self):
        with self.session("ui.choose('Pick', [('yes','Continue')])") as (process, _, __):
            self.receive(b'Continue')
            self.send('\x1b[200~q\nb\n2\n\x1b[201~')
            self.receive(timeout=.15)
            self.assertIsNone(process.poll())
            self.assertNotIn(b'RESULT=', self.received)
            self.send('\r')
            self.result('yes', process)

    def test_path_picker_menu_to_paste_preserves_quoted_path(self):
        with self.session("str(ui.pick_path('Choose DNA'))") as (process, _, __):
            self.receive(b'Paste or drag a path instead')
            self.send('\x1b[B\r')
            self.receive(b'Path (b to go back): ')
            self.send("\x1b[200~'/tmp/invented sample.txt'\n\x1b[201~")
            self.receive(timeout=.15)
            self.assertIsNone(process.poll())
            self.assertNotIn(b'RESULT=', self.received)
            self.send('\r')
            self.result(str(Path('/tmp/invented sample.txt').resolve()), process)

    def test_fresh_step_in_real_terminal_preserves_previous_output(self):
        expression="(ui.say('Previous step'), ui.screen(), ui.line('Next step'))[-1]"
        with self.session(expression) as (process, _, __):
            self.receive(b'Next step: ')
            self.assertIn(b'Previous step',self.received)
            self.assertIn(b'\x1b[H',self.received)
            for erase in (b'\x1b[3J', b'\x1b[2J', b'\x1b[?1049h', b'\x1bc'):
                self.assertNotIn(erase,self.received)
            self.send('finished\r')
            self.result('finished',process)

    def test_cancel_during_analysis_stops_tool_and_cleans_work(self):
        import signal
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);dna=root/'invented.txt';dna.write_bytes(DATA)
            pid_file=root/'tool.pid';out=root/'reports'
            child_code=f"import os,time;from pathlib import Path;Path({str(pid_file)!r}).write_text(str(os.getpid()));time.sleep(60)"
            setup=("import sys,tempfile\nfrom pathlib import Path\nfrom unittest.mock import patch\n"
                   "from genetics_app.__main__ import main\nfrom genetics_app.extras import run_command\n"
                   "def fake_pipeline(*args,**kwargs):\n"
                   f"    with tempfile.TemporaryDirectory(prefix='calculation-',dir={str(root)!r}) as work:\n"
                   f"        return run_command([sys.executable,'-c',{child_code!r}],work)\n"
                   "patch('genetics_app.terminal.reference_options',return_value=[('full','Full')]).start()\n"
                   "patch('genetics_app.pipeline.full_report',side_effect=fake_pipeline).start()\n")
            expression=f"main(['run',{str(dna)!r},'--refs',{str(root)!r},'--out',{str(out)!r}])"
            with self.session(expression,setup) as (process,slave,previous):
                self.receive(b'Save reports in')
                self.send('\r')
                self.receive(b'Review and generate report')
                self.send('\r')
                deadline=time.monotonic()+5
                while not pid_file.exists() and time.monotonic()<deadline:self.receive(timeout=.05)
                self.assertTrue(pid_file.exists())
                tool_pid=int(pid_file.read_text())
                os.kill(process.pid,signal.SIGINT)
                self.result(130,process)
                self.assert_restored(slave,previous)
                with self.assertRaises(ProcessLookupError):os.kill(tool_pid,0)
                self.assertFalse(list(root.glob('calculation-*')))
                self.assertFalse(list(out.glob('report.*')))

    def test_input_cancel_restores_terminal_even_inside_unfinished_paste(self):
        import signal
        with self.session() as (process, slave, previous):
            self.receive(b'Path: ')
            self.send('\x1b[200~/tmp/unfinished')
            os.kill(process.pid, signal.SIGINT)
            self.receive(b'KeyboardInterrupt')
            process.wait(timeout=3)
            self.receive(timeout=.1)
            self.assert_restored(slave, previous)
            self.assertIn(b'\x1b[?2004l', self.received)

    def test_long_paste_stays_editable_and_home_end_work(self):
        path = '/tmp/' + 'sample/'*30 + 'dna.txt'
        with self.session() as (process, _, __):
            self.receive(b'Path: ')
            self.send('\x1b[200~'+path+'\x1b[201~')
            self.send('\x01X\x04\x05\r')
            self.result('X'+path[1:], process)

    def test_defaults_and_redirected_input_are_preserved(self):
        with self.session("ui.line('Path', '/tmp/suggested')") as (process, _, __):
            self.receive(b'\x1b[?2004h')
            self.send('\r')
            self.result('/tmp/suggested', process)
        ui = Terminal(io.StringIO('/tmp/fixture\n'), io.StringIO())
        self.assertEqual(ui.line('Path'), '/tmp/fixture')


class BrowserMenuTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();self.root = Path(self.temp.name)
        self.output = io.StringIO();self.ui = Terminal(io.StringIO(), self.output)
        self.ui.arrows = True
        self.ui.line = lambda *args: self.fail('Browsing must not require typing a path')

    def tearDown(self):self.temp.cleanup()

    def test_browse_file_with_spaces_without_path_input(self):
        folder = self.root/'Samples';folder.mkdir()
        file = folder/'synthetic input.txt';file.write_bytes(DATA)
        with patch.object(self.ui, 'choose', side_effect=['browse', ('directory', folder), ('file', file)]):
            self.assertEqual(self.ui.pick_path('DNA', self.root/'missing.txt', suffixes={'.txt'}), file.resolve())

    def test_folder_navigation_and_parent(self):
        child = self.root/'child';child.mkdir()
        with patch.object(self.ui, 'choose', side_effect=['parent', 'use']):
            self.assertEqual(self.ui.browse(child, directory=True), self.root.resolve())

    def test_back_in_nested_folder_returns_to_the_previous_folder(self):
        child=self.root/'child';child.mkdir()
        with patch.object(self.ui,'choose',side_effect=[('directory',child),'back','use']) as choose:
            self.assertEqual(self.ui.browse(self.root,directory=True),self.root.resolve())
        self.assertEqual([call.args[0] for call in choose.call_args_list],['Folder: '+str(self.root.resolve()),'Folder: '+str(child.resolve()),'Folder: '+str(self.root.resolve())])

    def test_reexport_back_from_folder_returns_to_report_picker(self):
        from genetics_app.terminal import saved_report_flow
        source=export_report(report(),self.root/'saved',('json',))['json']
        with patch.object(self.ui,'pick_path',side_effect=[source,None,None]) as pick:
            self.assertIsNone(saved_report_flow(self.ui))
        self.assertEqual([call.args[0] for call in pick.call_args_list],['Choose a saved report','Export folder','Choose a saved report'])

    def test_invalid_saved_report_retries_in_the_file_picker(self):
        from genetics_app.terminal import saved_report_flow
        bad=self.root/'invalid.json';bad.write_text('{broken')
        source=export_report(report(),self.root/'saved',('json',))['json']
        with patch.object(self.ui,'pick_path',side_effect=[bad,source,self.root/'new']) as pick,patch.object(self.ui,'choose',return_value='done'):
            self.assertEqual(saved_report_flow(self.ui),0)
        self.assertIn('Could not open this report',pick.call_args_list[1].kwargs['notice'])

    def test_hidden_files_toggle_and_file_type_filter(self):
        hidden = self.root/'.synthetic.txt';hidden.write_bytes(DATA)
        (self.root/'unrelated.json').write_text('{}')
        def choose(prompt, options, **kwargs):
            labels = [label for key, label in options]
            if 'Show hidden items' in labels:
                self.assertNotIn(hidden.name, labels);self.assertNotIn('unrelated.json', labels)
                return 'hidden'
            self.assertIn(hidden.name, labels)
            return ('file', hidden)
        with patch.object(self.ui, 'choose', side_effect=choose):
            self.assertEqual(self.ui.browse(self.root, suffixes={'.txt'}), hidden.resolve())

    def test_suggested_output_needs_no_typing_or_premature_creation(self):
        suggested = self.root/'new report folder'
        with patch.object(self.ui, 'choose', return_value='suggested'):
            self.assertEqual(self.ui.pick_path('Save', suggested, directory=True, create=True), suggested.resolve())
        self.assertFalse(suggested.exists())

    def test_manual_path_starts_blank_and_empty_enter_does_not_select_a_folder(self):
        suggested = self.root/'suggested reports'
        chosen = self.root/'chosen reports'
        with patch.object(self.ui, 'choose', return_value='paste'), patch.object(self.ui, 'line', side_effect=['', str(chosen)]) as line:
            self.assertEqual(self.ui.pick_path('Save', suggested, directory=True, create=True), chosen.resolve())
        self.assertEqual(line.call_count, 2)
        for call in line.call_args_list:
            self.assertEqual(call.args, ('Path (b to go back)',))
            self.assertEqual(call.kwargs, {})
        self.assertIn('Enter a path', self.output.getvalue())
        self.assertFalse(chosen.exists())

    def test_browse_back_returns_to_picker_then_back_to_previous_step(self):
        with patch.object(self.ui, 'choose', side_effect=['browse', 'back', 'back']):
            self.assertIsNone(self.ui.pick_path('DNA', self.root/'missing.txt'))

    def test_long_directory_menu_is_scrolled_instead_of_dumped(self):
        options = [(str(i), 'File '+str(i)) for i in range(50)]
        with patch.object(self.ui, '_key', side_effect=['down']*15 + ['\n']):
            self.assertEqual(self.ui.choose_keys(options, 0, True), '15')
        text = self.output.getvalue()
        self.assertIn('of 50', text);self.assertNotIn('File 49', text)

    def test_home_export_opens_existing_json_without_analysis_or_commands(self):
        from genetics_app.terminal import application
        source = self.root/'existing';source.mkdir()
        paths = export_report(report(), source, ('json',))
        out = self.root/'exports'
        with patch.object(self.ui,'choose',side_effect=['export','all','done']), patch.object(self.ui,'pick_path',side_effect=[paths['json'],out]), patch('genetics_app.pipeline.full_report',side_effect=AssertionError('No rerun')):
            self.assertEqual(application(terminal=self.ui),0)
        self.assertTrue((out/'report.md').exists());self.assertTrue((out/'report.txt').exists())
        self.assertEqual(json.loads((out/'report.json').read_text()),json.loads(paths['json'].read_text()))

    def test_home_example_all_formats_without_typing(self):
        from genetics_app.terminal import application
        out = self.root/'example'
        with patch.object(self.ui,'choose',side_effect=['example','all','done']), patch.object(self.ui,'pick_path',return_value=out):
            self.assertEqual(application(refs=self.root/'no-refs',terminal=self.ui),0)
        self.assertTrue(json.loads((out/'report.json').read_text())['qc']['synthetic_input_declared'])
        self.assertTrue((out/'report.md').exists());self.assertTrue((out/'report.txt').exists())

    def test_recovery_menu_does_not_offer_reduced_dna_analysis(self):
        from genetics_app.terminal import application, QuitFlow
        with patch.object(self.ui,'choose',side_effect=QuitFlow) as choose:
            with self.assertRaises(QuitFlow):application(terminal=self.ui,limited=True)
        options=choose.call_args.args[1]
        self.assertEqual([key for key,label in options],['export','example'])
        self.assertIn('Full setup is not being used',self.output.getvalue())

    def test_recovery_example_ignores_full_profile_and_uses_no_research_packages(self):
        from genetics_app.terminal import application
        (self.root/'pipeline.json').write_text('{invalid full profile fixture}')
        out=self.root/'recovery-example'
        with patch.object(self.ui,'choose',side_effect=['example','all','done']), patch.object(self.ui,'pick_path',return_value=out), patch('genetics_app.pipeline.full_report',side_effect=AssertionError('No research dependencies')):
            self.assertEqual(application(refs=self.root,terminal=self.ui,limited=True),0)
        self.assertEqual(len(json.loads((out/'report.json').read_text())['traits']),3)

    def test_back_from_initial_picker_returns_to_home(self):
        from genetics_app.terminal import application, QuitFlow
        with patch.object(self.ui,'choose',side_effect=['analyze',QuitFlow]), patch.object(self.ui,'pick_path',return_value=None):
            with self.assertRaises(QuitFlow):application(terminal=self.ui)


class LetterFormatTests(unittest.TestCase):
    def test_common_explicit_headers(self):
        from genetics_app.parser import detect_strand
        cases = [
            ('# genotype call oriented with respect to the plus strand on the human reference sequence.', 'forward'),
            ('# Genotypes are reported\n# on the forward (+) strand with respect to the human reference.', 'forward'),
            ('# GRCh37 forward-strand calls at published marker coordinates.', 'forward'),
            ('# strand: +', 'forward'), ('# orientation=minus', 'reverse'),
            ('# Data reported on the reverse (-) strand.', 'reverse')]
        for header, expected in cases:
            with self.subTest(header=header):
                direction, evidence = detect_strand((header+'\nrs1 1 100 AG').encode())
                self.assertEqual(direction, expected);self.assertTrue(evidence)

    def test_unknown_and_unrelated_prose_do_not_guess(self):
        from genetics_app.parser import detect_strand
        for header in ['# 23andMe build 37', '# Gene ABC is on the minus strand',
                       '# Data may be on the forward strand', '# Data is not on the forward strand',
                       '# Data: forward strand not confirmed', '# build 37']:
            with self.subTest(header=header):
                self.assertEqual(detect_strand((header+'\nrs1 1 100 AG').encode()), ('unknown', []))

    def test_negated_reverse_does_not_cancel_explicit_forward(self):
        from genetics_app.parser import detect_strand
        self.assertEqual(detect_strand(b'# Genotypes are reported on the forward strand, not the reverse strand.')[0], 'forward')

    def test_conflicting_direction_metadata_rejected(self):
        from genetics_app.parser import detect_strand
        with self.assertRaisesRegex(ValueError, 'both forward/plus and reverse/minus'):
            detect_strand(b'# strand: plus\n# strand: minus')

    def test_parser_rejects_reverse_even_with_forward_argument(self):
        from genetics_app.parser import parse_genome
        with self.assertRaisesRegex(ValueError, 'reverse/minus'):
            parse_genome(b'# build 37\n# Data is reported on the reverse strand.\nrs1 1 100 AG', strand='forward')

    def test_detected_format_needs_no_prompt_or_user_action(self):
        from genetics_app.terminal import confirm_letter_format
        output=io.StringIO();ui=Terminal(io.StringIO('\n'),output)
        with patch.object(ui,'choose',side_effect=AssertionError('Detected format needs no confirmation')):
            self.assertEqual(confirm_letter_format(DATA,ui),'yes')
        text=output.getvalue()
        self.assertEqual(text, '')

    def test_unknown_defaults_to_help_instead_of_confirmation(self):
        from genetics_app.terminal import confirm_letter_format
        output=io.StringIO();ui=Terminal(io.StringIO('\n'),output)
        self.assertEqual(confirm_letter_format(b'# build 37\nrs1 1 100 AG',ui),'help')
        self.assertIn("I'm not sure - show me how to check",output.getvalue())

    def test_unknown_help_can_return_to_file_without_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);dna=root/'invented.txt';dna.write_bytes(DATA.replace(b'# Genotype calls use the forward strand.\n',b''))
            ui=Terminal(io.StringIO(str(dna)+'\n\n\nb\n'),io.StringIO())
            with patch('genetics_app.terminal.reference_options',return_value=[('full','Full analysis')]),patch('genetics_app.pipeline.full_report') as run:
                self.assertIsNone(guided(root,root/'reports',terminal=ui))
            run.assert_not_called();self.assertFalse((root/'reports/report.json').exists())
            self.assertIn('How to check',ui.stdout.getvalue())

    def test_short_menu_spacing_and_no_meaningless_page_counter(self):
        output=io.StringIO();ui=Terminal(io.StringIO(),output);ui.arrows=True
        with patch.object(ui,'key_mode'),patch.object(ui,'_key',side_effect=['down','\n']):
            self.assertEqual(ui.choose('DNA letter format',[('yes','Continue'),('help','Explain')]),'help')
        text=output.getvalue()
        self.assertIn('DNA letter format\n'+'─'*52+'\n\n',text)
        self.assertIn('↑/↓ move   Enter select\n\n',text);self.assertNotIn('of 4',text)
        self.assertNotIn('b back',text);self.assertNotIn('q quit',text)
        self.assertIn('Back',text);self.assertIn('Quit',text)
        self.assertIn('\x1b[4A',text)

    def test_narrow_menu_labels_fit_the_terminal(self):
        import os
        output=io.StringIO();ui=Terminal(io.StringIO(),output)
        with patch('genetics_app.terminal.shutil.get_terminal_size',return_value=os.terminal_size((20,16))),patch.object(ui,'_key',return_value='\n'):
            self.assertEqual(ui.choose_keys([('yes','A long label with spaces and a long ending')],0,False),'yes')
        lines=output.getvalue().splitlines()
        rendered=[line.replace('\x1b[2K','') for line in lines if '> ' in line]
        self.assertTrue(rendered);self.assertTrue(all(len(line)<=20 for line in rendered))
