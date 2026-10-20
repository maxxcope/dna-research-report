"""Guided terminal input, back navigation and export selection."""
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
import os
import shlex
import shutil
import sqlite3
import sys
from .local import DEFAULT_REFS, read_input, analyze_local
from .parser import detect_build, detect_strand, parse_genome
from .progress import Progress
from .report import export_report, safe_text, output_folder, read_report


class QuitFlow(Exception):
    pass


def terminal_size():
    """Python 3.10 can return zero dimensions for a newly opened PTY."""
    size = shutil.get_terminal_size((80, 24))
    return os.terminal_size((size.columns if size.columns > 0 else 80,
                             size.lines if size.lines > 0 else 24))


class Terminal:
    def __init__(self, stdin=None, stdout=None):
        self.stdin = stdin or sys.stdin
        self.stdout = stdout or sys.stdout
        self.arrows = self.stdin.isatty() and self.stdout.isatty()

    def say(self, text=''):
        self.stdout.write(safe_text(text) + '\n')
        self.stdout.flush()

    def screen(self):
        """Start a fresh page by scrolling, never by erasing terminal history."""
        if self.arrows:
            rows = max(1, terminal_size().lines)
            self.stdout.write('\r\n' * rows + '\x1b[H')
            self.stdout.flush()

    def line(self, prompt, default=''):
        suffix = f' [{safe_text(default)}]' if default else ''
        label = safe_text(prompt) + suffix + ': '
        if self.arrows and os.name != 'nt':
            value = self.edit_line(label)
        else:
            self.stdout.write(label)
            self.stdout.flush()
            value = self.stdin.readline()
            if not value:
                raise QuitFlow()
        value = value.strip()
        if value.lower() == 'q':
            raise QuitFlow()
        return value or str(default)

    def edit_line(self, label):
        """Keep terminal pastes editable; only a separate Enter submits them."""
        text = ''
        cursor = 0
        with self.key_mode():
            while True:
                # Keep long paths on one screen row so editing never corrupts
                # the next prompt. The full value remains in the input buffer.
                width = max(4, terminal_size().columns)
                prompt = label if len(label) < width // 2 else label[:max(0, width // 2 - 3)] + '…: '
                room = max(1, width - len(prompt) - 1)
                start = max(0, cursor - room + 1)
                visible = safe_text(text[start:start + room])
                self.stdout.write('\r\x1b[2K' + prompt + visible)
                retreat = len(visible) - (cursor - start)
                if retreat > 0:
                    self.stdout.write(f'\x1b[{retreat}D')
                self.stdout.flush()
                key = self._key()
                if isinstance(key, tuple):
                    pasted = key[1].replace('\r\n', '\n').replace('\r', '\n').strip('\n')
                    # Newlines and terminal controls are data, never commands.
                    pasted = ''.join(' ' if char in '\n\t' else char for char in pasted
                                     if char.isprintable() or char in '\n\t')
                    text = text[:cursor] + pasted + text[cursor:]
                    cursor += len(pasted)
                elif key in ('\r', '\n'):
                    self.stdout.write('\n')
                    return text
                elif key == '\x03':
                    raise KeyboardInterrupt()
                elif key == '\x04' and not text:
                    raise QuitFlow()
                elif key in ('\x7f', '\b') and cursor:
                    text = text[:cursor - 1] + text[cursor:]
                    cursor -= 1
                elif key in ('delete', '\x04'):
                    text = text[:cursor] + text[cursor + 1:]
                elif key in ('left', '\x02'):
                    cursor = max(0, cursor - 1)
                elif key in ('right', '\x06'):
                    cursor = min(len(text), cursor + 1)
                elif key in ('home', '\x01'):
                    cursor = 0
                elif key in ('end', '\x05'):
                    cursor = len(text)
                elif key == '\x15':
                    text, cursor = text[cursor:], 0
                elif key == '\x0b':
                    text = text[:cursor]
                elif len(key) == 1 and key.isprintable():
                    text = text[:cursor] + key + text[cursor:]
                    cursor += 1

    def _character(self):
        import codecs
        decoder = codecs.getincrementaldecoder('utf-8')(errors='replace')
        while True:
            byte = os.read(self.stdin.fileno(), 1)
            if not byte:
                raise QuitFlow()
            char = decoder.decode(byte)
            if char:
                return char

    def _paste(self):
        data = bytearray()
        ending = b'\x1b[201~'
        while not data.endswith(ending):
            byte = os.read(self.stdin.fileno(), 1)
            if not byte:
                raise QuitFlow()
            data.extend(byte)
        return ('paste', data[:-len(ending)].decode('utf-8', errors='replace'))

    def _key(self):
        if os.name == 'nt':
            import msvcrt
            key = msvcrt.getwch()
            if key in ('\x00', '\xe0'):
                return {'H': 'up', 'P': 'down'}.get(msvcrt.getwch(), '')
            return key
        import select
        fd = self.stdin.fileno()
        key = self._character()
        if key == '\x1b' and select.select([fd], [], [], .08)[0]:
            tail = os.read(fd, 1)
            if tail not in (b'[', b'O'):
                return ''
            while select.select([fd], [], [], .08)[0]:
                byte = os.read(fd, 1)
                if not byte:
                    raise QuitFlow()
                tail += byte
                if 0x40 <= byte[0] <= 0x7e:
                    if tail == b'[200~':
                        return self._paste()
                    return {b'[A': 'up', b'OA': 'up', b'[B': 'down', b'OB': 'down',
                            b'[C': 'right', b'OC': 'right', b'[D': 'left', b'OD': 'left',
                            b'[H': 'home', b'OH': 'home', b'[1~': 'home', b'[7~': 'home',
                            b'[F': 'end', b'OF': 'end', b'[4~': 'end', b'[8~': 'end',
                            b'[3~': 'delete'}.get(tail, '')
                if len(tail) > 32:
                    return ''
            return ''
        return key

    @contextmanager
    def key_mode(self):
        if os.name == 'nt':
            yield
            return
        import termios
        import tty
        fd = self.stdin.fileno()
        previous = termios.tcgetattr(fd)
        try:
            # Do not flush fast keystrokes between menu redraws.
            tty.setcbreak(fd, termios.TCSANOW)
            current = termios.tcgetattr(fd)
            # Preserve pasted CRLF pairs so they normalize to one separator.
            current[0] &= ~termios.ICRNL
            termios.tcsetattr(fd, termios.TCSANOW, current)
            self.stdout.write('\x1b[?2004h')
            self.stdout.flush()
            yield
        finally:
            try:
                self.stdout.write('\x1b[?2004l')
                self.stdout.flush()
            finally:
                termios.tcsetattr(fd, termios.TCSADRAIN, previous)

    def choose(self, prompt, options, default=0, back=True):
        options = list(options)
        if back:
            options.append(('back', 'Back'))
        options.append(('quit', 'Quit'))
        self.say('\n' + prompt)
        self.say('─' * min(52, max(12, terminal_size().columns - 2)))
        self.say()
        if not self.arrows:
            for i, (_, label) in enumerate(options, 1):
                self.say(f'  {i}. {label}')
            while True:
                choice = self.line('Choose a number', str(default + 1))
                if choice.lower() == 'b' and back:
                    return 'back'
                if choice.isdigit() and 1 <= int(choice) <= len(options):
                    result = options[int(choice) - 1][0]
                    if result == 'quit':
                        raise QuitFlow()
                    return result
                self.say('Choose one of the listed numbers.')
        with self.key_mode():
            return self.choose_keys(options, default, back)

    def choose_keys(self, options, default, back):
        current = default
        self.say('↑/↓ move   Enter select')
        self.say()
        size = terminal_size()
        rows = min(len(options), max(1, min(12, size.lines - 9)))
        columns = max(1, size.columns - 4)
        scrolling = len(options) > rows
        while True:
            start = min(max(0, current - rows // 2), len(options) - rows)
            if scrolling:
                self.stdout.write('\r\x1b[2K')
                self.say(f'{start + 1}-{start + rows} of {len(options)}')
            for i in range(start, start + rows):
                label = safe_text(options[i][1]).replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
                if len(label) > columns:
                    if columns >= 8:
                        head = columns // 3
                        label = label[:head] + '...' + label[-(columns - head - 3):]
                    else:
                        label = label[:columns]
                self.stdout.write('\r\x1b[2K')
                self.say(('  > ' if i == current else '    ') + label)
            key = self._key()
            if isinstance(key, tuple):
                # A paste over a selector must not trigger shortcuts or Enter.
                key = ''
            if key in ('\x03',):
                raise KeyboardInterrupt()
            if key.lower() == 'q' or key == '\x04':
                raise QuitFlow()
            if key in ('b', 'B', '\x1b') and back:
                return 'back'
            if key in ('\r', '\n'):
                result = options[current][0]
                if result == 'quit':
                    raise QuitFlow()
                return result
            if key in ('up', 'k'):
                current = (current - 1) % len(options)
            elif key in ('down', 'j'):
                current = (current + 1) % len(options)
            elif key.isdigit() and 1 <= int(key) <= len(options):
                current = int(key) - 1
            self.stdout.write(f'\x1b[{rows + int(scrolling)}A\r')
            self.stdout.flush()

    def pick_path(self, prompt, suggested=None, directory=False, suffixes=None, create=False, notice=None):
        """Browse by default in a terminal; retain line input for redirected flows."""
        if not self.arrows:
            if notice:
                self.say(notice)
            text = self.line(prompt + ' (b to go back)', str(suggested) if suggested else '')
            return None if text.lower() == 'b' else path_value(text)
        suggested = Path(suggested).expanduser().resolve() if suggested else None
        options = []
        if suggested and ((directory and (suggested.is_dir() or create)) or (not directory and suggested.is_file())):
            options.append(('suggested', 'Use ' + suggested.name + (' (suggested folder)' if directory else '')))
        options += [('browse', 'Browse folders' if directory else 'Browse files'), ('paste', 'Paste or drag a path instead')]
        while True:
            self.screen()
            if notice:
                self.say(notice)
                notice = None
            if options[0][0] == 'suggested':
                self.say('Suggested location: ' + str(suggested))
            choice = self.choose(prompt, options)
            if choice == 'back':
                return None
            if choice == 'suggested':
                return suggested
            if choice == 'paste':
                self.screen()
                self.say('Paste or type the path, edit it if needed, then press Enter.')
                while True:
                    text = self.line('Path (b to go back)')
                    if text.lower() == 'b':
                        break
                    if text.strip():
                        return path_value(text)
                    self.say('Enter a path, or type b to go back.')
            else:
                start = suggested if directory else suggested.parent if suggested else Path.home()
                while start and not start.is_dir():
                    start = start.parent
                result = self.browse(start or Path.home(), directory, suffixes, create)
                if result is not None:
                    return result

    def browse(self, start, directory=False, suffixes=None, create=False):
        current = Path(start).expanduser().resolve()
        hidden = False
        history = []
        while True:
            self.screen()
            options = []
            if directory:
                options.append(('use', 'Use this folder'))
                if create:
                    options.append(('new', 'Create a dated report folder here'))
            if current.parent != current:
                options.append(('parent', '../  Parent folder'))
            places = [('Home', Path.home()), ('Downloads', Path.home()/'Downloads'),
                      ('Desktop', Path.home()/'Desktop'), ('Documents', Path.home()/'Documents')]
            for label, path in places:
                if path.is_dir() and path.resolve() != current:
                    options.append((('place', path), 'Go to ' + label))
            options.append(('hidden', 'Hide hidden items' if hidden else 'Show hidden items'))
            try:
                entries = sorted(current.iterdir(), key=lambda p: (not p.is_dir(), p.name.casefold()))
                for path in entries:
                    if not hidden and path.name.startswith('.'):
                        continue
                    if path.is_dir():
                        options.append((('directory', path), path.name + '/'))
                    elif not directory and (not suffixes or path.suffix.lower() in suffixes):
                        options.append((('file', path), path.name))
            except OSError as error:
                self.say('Cannot read this folder: ' + str(error))
            choice = self.choose(('Folder: ' if directory else 'Files: ') + str(current), options)
            if choice == 'back':
                if history:
                    current = history.pop()
                    continue
                return None
            if choice == 'use':
                return current
            if choice == 'new':
                candidate = current / default_output().name
                index = 2
                while candidate.exists():
                    candidate = current / (default_output().name + '-' + str(index))
                    index += 1
                return candidate
            if choice == 'parent':
                history.append(current)
                current = current.parent
            elif choice == 'hidden':
                hidden = not hidden
            elif choice[0] == 'file':
                return choice[1].resolve()
            else:
                history.append(current)
                current = choice[1].resolve()


def path_value(text):
    value = text.strip()
    if not value:
        raise ValueError('Enter a file or folder path.')
    path = Path(value).expanduser()
    if path.exists():
        return path.resolve()
    try:
        parts = shlex.split(value)
    except ValueError:
        parts = []
    if len(parts) == 1:
        if not parts[0].strip():
            raise ValueError('Enter a file or folder path.')
        return Path(parts[0]).expanduser().resolve()
    return path.resolve()


def default_output():
    return Path.home() / 'dna-reports' / datetime.now().strftime('%Y%m%d-%H%M%S')


def reference_options(refs, build):
    """Inspect public metadata only; never silently reduce the analysis."""
    refs = Path(refs)
    options = []
    profile = refs / 'pipeline.json'
    if profile.is_file():
        from .pipeline import load_config
        from .references import metadata
        config = load_config(profile)
        if metadata(config['clinvar']).get('build') == build:
            options.append(('full', 'Full analysis (recommended): findings, audits and additional analyses'))
    if (refs / f'clinvar{build}.sqlite').is_file():
        options.append(('basic', 'Limited analysis: clinical annotations and three curated traits'))
    return options


def describe_result(report, terminal):
    pipeline = report.get('pipeline')
    if pipeline and pipeline['status'] == 'completed':
        terminal.say('\nSuccess. All configured analysis stages completed.')
    elif pipeline:
        terminal.say('\nAnalysis finished with incomplete or unavailable stages:')
        for stage in pipeline['stages']:
            if stage['status'] != 'completed':
                terminal.say('  ' + stage['name'] + ': ' + stage['status'].replace('_', ' '))
    else:
        terminal.say('\nLimited report complete. The full research pipeline was not run.')
    counts = report.get('knowledge', {}).get('counts', {})
    if counts:
        terminal.say(f"{counts['variant_summaries']:,} genotype explanations; {counts['genosets']} multi-marker findings; {counts['unresolved_loci']:,} unresolved source loci.")
    terminal.say(f"{len(report['clinvar']):,} clinical annotations; {len(report['gwas']):,} study associations.")


def export_menu(report, out, terminal, saved, notice=None):
    written = dict(saved)
    all_saved = list(saved.items())
    exported = False
    try:
        while True:
            terminal.screen()
            describe_result(report, terminal)
            if notice:
                terminal.say(notice)
                notice = None
            terminal.say('\nSaved files:')
            for path in written.values():
                terminal.say('  ' + str(path))
            if 'bundle' in written:
                terminal.say('Start with: ' + str(written['bundle'] / 'overview.md'))
            labels = {'bundle': 'Organized report folder with overview, parts and full archive (recommended)',
                      'json': 'Full uncompressed JSON (.json)',
                      'markdown': 'Full single Markdown file (.md)', 'text': 'Full single text file (.txt)'}
            options = [(fmt, label) for fmt, label in labels.items() if fmt not in written]
            if len(options) > 1:
                options.append(('all', 'All formats'))
            options += [('folder', 'Export to another folder'), ('done', 'Finish')]
            preferred = 'done' if exported or 'bundle' in written else 'bundle'
            choice = terminal.choose('Report saved. Export another format?' if 'bundle' in written else 'Export this completed analysis', options,
                                     default=next(i for i, (key, _) in enumerate(options) if key == preferred), back=False)
            if choice == 'done':
                return
            if choice == 'folder':
                folder = terminal.pick_path('Choose export folder', out, directory=True, create=True)
                if folder is not None:
                    out = folder
                    written = {fmt: path for fmt, path in all_saved
                               if path.parent.resolve() == Path(folder).resolve()}
                    exported = False
                continue
            formats = tuple(labels) if choice == 'all' else (choice,)
            formats = tuple(f for f in formats if f not in written)
            if not formats:
                notice = 'That export is already saved. Choose another format or folder.'
                continue
            try:
                terminal.say('Writing complete export...')
                new_files = export_report(report, out, formats)
                written.update(new_files)
                all_saved.extend(new_files.items())
                exported = True
            except (OSError, ValueError) as error:
                notice = 'Export failed: ' + str(error) + '\nThe completed analysis is still available. Choose another export folder to retry.'

    except (QuitFlow, KeyboardInterrupt) as error:
        if all_saved:
            terminal.say('\nSaved reports:')
            for fmt, path in dict.fromkeys(all_saved):
                terminal.say('  ' + str(path))
        else:
            terminal.say('No report was saved.')
        if isinstance(error, KeyboardInterrupt):
            raise


def confirm_letter_format(data, terminal):
    """Explain orientation using file evidence and a useful unknown path."""
    direction, _ = detect_strand(data)
    if direction == 'forward':
        return 'yes'
    terminal.say('\nDNA letters can be written in two directions. This app uses the format called "forward" or "plus" to match them correctly.')
    if direction == 'reverse':
        terminal.say('This file uses the opposite format (reverse/minus). Ask the provider for a forward/plus export; editing the header will not convert it.')
        return terminal.choose('Choose a different export', [('file', 'Choose another DNA file')])
    terminal.say('\nNo direction label was found in the file. Check before continuing so the DNA letters are interpreted correctly.')
    return terminal.choose('Check your file format', [
        ('help', "I'm not sure - show me how to check"),
        ('yes', 'I checked: this file uses forward/plus'),
        ('file', 'Choose another DNA file')])


def explain_letter_format(terminal):
    terminal.say('\nHow to check')
    terminal.say('1. Open the DNA text file in a text editor. Read the comment lines at the top, usually starting with #.')
    terminal.say('2. Look for "forward strand" or "plus strand". Both mean the format this app needs.')
    terminal.say('3. If neither is shown, check the provider\'s raw-data download instructions or ask their support team. Do not guess.')
    terminal.say('\nA build number such as GRCh37 or GRCh38 describes the DNA map. It does not tell you which direction the letters use.')


def select_references(refs, build, terminal):
    """Use a matching full profile automatically; never silently use less."""
    interacted = False
    notice = None
    while True:
        try:
            options = reference_options(refs, build)
        except (OSError, ValueError, KeyError, TypeError, AttributeError, sqlite3.Error) as error:
            options = []
            notice = 'Could not read these references: ' + str(error)
        if any(key == 'full' for key, _ in options):
            return Path(refs), 'full', interacted
        interacted = True
        terminal.screen()
        if notice:
            terminal.say(notice)
            notice = None
        terminal.say(f'No full analysis database for GRCh{build} was found in:\n{refs}')
        terminal.say('Choose another prepared reference folder. The automatic installation prepares the full GRCh37 database.')
        choice = terminal.choose('Available analysis', options + [('folder', 'Choose another reference folder')])
        if choice == 'back':
            return None
        if choice == 'folder':
            folder = terminal.pick_path('Choose prepared reference folder', refs, directory=True)
            if folder is not None:
                refs = folder
        else:
            return Path(refs), choice, interacted


def guided(refs=DEFAULT_REFS, out=None, input_path=None, terminal=None):
    terminal = terminal or Terminal()
    state = 0
    history = []
    selected_path = str(input_path) if input_path else ''
    refs = Path(refs).expanduser()
    out = Path(out).expanduser() if out else default_output()
    notice = None
    editing_output = False
    editing_input = False
    build = 'auto'
    if input_path is not None:
        try:
            selected_path = str(Path(input_path).expanduser().resolve())
            data = read_input(selected_path)
            state = 1
            history = [0]
        except (OSError, ValueError) as error:
            notice = 'Cannot read this DNA file: ' + str(error)

    def back():
        return history.pop() if history else 0

    def remember(page):
        if not history or history[-1] != page:
            history.append(page)

    while True:
        try:
            if state not in (0, 2, 4):
                terminal.screen()
                terminal.say('DNA Research Report')
                if notice:
                    terminal.say(notice)
                    notice = None
            if state == 0:
                path = terminal.pick_path('Choose your DNA file', selected_path or None,
                                          suffixes={'.txt', '.csv', '.tsv'}, notice=notice)
                notice = None
                if path is None:
                    if editing_input:
                        state = 5
                        editing_input = False
                        continue
                    return None
                data = read_input(path)
                editing_input = False
                selected_path = str(path)
                history = [0]
                build = 'auto'
                state = 1
            elif state == 1:
                try:
                    detected = detect_build(data)
                    direction, _ = detect_strand(data)
                    if direction == 'reverse':
                        raise ValueError('This file uses reverse/minus DNA letters. Get a forward/plus export from the provider; changing the header does not convert it.')
                    if detected:
                        build = 'auto'
                    else:
                        terminal.say('Your file does not name its genome build: the version of the DNA map used for marker positions.')
                        terminal.say('Check the provider\'s raw-data download instructions for GRCh37/hg19 or GRCh38/hg38. Choosing the wrong map can match the wrong positions.')
                        choice = terminal.choose('Which map does your provider specify?', [
                            ('37', 'GRCh37 / hg19'), ('38', 'GRCh38 / hg38'),
                            ('back', 'I do not know - choose another file')], default=2, back=False)
                        if choice == 'back':
                            state = back()
                            continue
                        build = choice
                        remember(1)
                    genome = parse_genome(data, build)
                except ValueError as error:
                    notice = str(error)
                    state = 0
                    continue
                state = 2
            elif state == 2:
                if notice:
                    terminal.screen()
                    terminal.say(notice)
                    notice = None
                result = select_references(refs, genome.build, terminal)
                if result is None:
                    state = back()
                    continue
                refs, mode, interacted = result
                if interacted:
                    remember(2)
                state = 4 if direction == 'forward' else 3
            elif state == 3:
                choice = confirm_letter_format(data, terminal)
                if choice == 'back':
                    state = back()
                elif choice == 'file':
                    state = 0
                elif choice == 'help':
                    terminal.screen()
                    explain_letter_format(terminal)
                    answer = terminal.choose('Next step', [
                        ('file', 'Choose another DNA file'),
                        ('yes', 'I checked: this file uses forward/plus')])
                    if answer == 'file':
                        state = 0
                    elif answer == 'yes':
                        remember(3)
                        state = 4
                elif choice == 'yes':
                    remember(3)
                    state = 4
            elif state == 4:
                folder = terminal.pick_path('Save reports in', out, directory=True, create=True, notice=notice)
                notice = None
                if folder is None:
                    state = 5 if editing_output else back()
                    editing_output = False
                    continue
                candidate = output_folder(folder)
                if any((candidate / ('report.' + ext)).exists() or (candidate / ('report.' + ext)).is_symlink()
                       for ext in ('json', 'md', 'txt')):
                    raise ValueError('This folder already contains a report. Choose another folder or create a dated subfolder.')
                out = candidate
                if not editing_output:
                    remember(4)
                editing_output = False
                state = 5
            else:
                terminal.say(f'\nDNA file: {selected_path}\nGenome map: GRCh{genome.build} '
                             + ('(read from file)' if detected else '(selected by you)'))
                terminal.say('DNA format: forward/plus ' + ('(read from file)' if direction == 'forward' else '(confirmed by you)'))
                terminal.say('Analysis: ' + ('Full analysis' if mode == 'full' else 'Limited clinical and trait analysis'))
                terminal.say(f'References: {refs}\nSave reports in: {out}')
                choice = terminal.choose('Review and generate report', [
                    ('run', 'Generate report'), ('file', 'Change DNA file'),
                    ('output', 'Change report folder'), ('refs', 'Change reference folder')])
                if choice == 'back':
                    state = back()
                    continue
                if choice == 'file':
                    editing_input = True
                    state = 0
                    continue
                if choice == 'output':
                    editing_output = True
                    state = 4
                    continue
                if choice == 'refs':
                    folder = terminal.pick_path('Choose prepared reference folder', refs, directory=True)
                    if folder is not None:
                        result = select_references(folder, genome.build, terminal)
                        if result is not None:
                            refs, mode, _ = result
                    continue
                try:
                    output_folder(out, create=True)
                except (OSError, ValueError) as error:
                    notice = 'Cannot save reports here: ' + str(error)
                    editing_output = True
                    state = 4
                    continue
                try:
                    terminal.screen()
                    terminal.say('Generating report\nDNA file: ' + selected_path + '\nSave reports in: ' + str(out))
                    with Progress() as progress:
                        if mode == 'full':
                            from .pipeline import full_report
                            report = full_report(data, refs / 'pipeline.json', build, events=progress.event)
                        else:
                            progress.event(dict(kind='detail', message='Matching clinical annotations and curated traits'))
                            report = analyze_local(data, build, refs)
                except (OSError, ValueError, sqlite3.Error) as error:
                    notice = 'Analysis failed: ' + str(error) + '\nNo report was exported. Change the file or references, or retry.'
                    continue
                export_notice = None
                try:
                    terminal.say('Saving report maps, findings and full archive...')
                    saved = export_report(report, out, ('bundle',))
                except (OSError, ValueError) as error:
                    export_notice = 'Could not save report folder: ' + str(error) + '\nChoose another export folder to save the completed analysis.'
                    saved = {}
                export_menu(report, out, terminal, saved, notice=export_notice)
                return 0
        except (OSError, ValueError) as error:
            notice = 'Cannot continue: ' + str(error)



def saved_report_flow(terminal):
    source = None
    report = None
    notice = None
    while True:
        if report is None:
            source = terminal.pick_path('Choose a saved report', source, suffixes={'.json', '.gz'}, notice=notice)
            notice = None
            if source is None:
                return None
            try:
                report = read_report(source)
            except (OSError, ValueError) as error:
                notice = 'Could not open this report: ' + str(error)
                continue
        base = source.parent.parent if source.name == 'full-report.json.gz' and source.parent.name == 'report-bundle' else source.parent
        folder = terminal.pick_path('Export folder', base, directory=True, create=True)
        if folder is None:
            report = None
            continue
        saved = {'json': source} if source.suffix == '.json' and source.parent.resolve() == Path(folder).resolve() else {}
        if source.name == 'full-report.json.gz' and source.parent.name == 'report-bundle' and source.parent.parent.resolve() == Path(folder).resolve():
            saved['bundle'] = source.parent
        export_menu(report, folder, terminal, saved)
        return 0


def application(refs=DEFAULT_REFS, terminal=None, limited=False):
    """Normal entry point: all everyday actions are selected from menus."""
    terminal = terminal or Terminal()
    notice = None
    while True:
        terminal.screen()
        terminal.say('DNA Research Report')
        if limited:
            terminal.say('Full setup is not being used. Export a saved report or try the three-rule example.')
        if notice:
            terminal.say(notice)
            notice = None
        options = [('export', 'Export a saved report'), ('example', 'Try an invented example')]
        if not limited:
            options.insert(0, ('analyze', 'Analyze a DNA file'))
        action = terminal.choose('What would you like to do?', options, back=False)
        if action == 'analyze':
            result = guided(refs, terminal=terminal)
            if result is not None:
                return result
            continue
        if action == 'export':
            result = saved_report_flow(terminal)
            if result is not None:
                return result
            continue
        else:
            from importlib.resources import files
            folder = terminal.pick_path('Save example reports in', default_output(), directory=True, create=True)
            if folder is None:
                continue
            try:
                data = files('genetics_app').joinpath('data/demo_dna.txt').read_bytes()
                terminal.screen()
                terminal.say('Generating example report')
                if not limited and (Path(refs)/'pipeline.json').is_file():
                    from .pipeline import full_report
                    with Progress() as progress:
                        report = full_report(data, Path(refs)/'pipeline.json', events=progress.event)
                else:
                    from .engine import analyze
                    terminal.say('This example demonstrates three curated trait rules; the full research pipeline is not being used.')
                    report = analyze(data, traits_only=True)
                saved = export_report(report, folder, ('bundle',))
                try:
                    export_menu(report, folder, terminal, saved)
                except QuitFlow:
                    terminal.say('Complete report folder saved: ' + str(saved['bundle']))
                return 0
            except (OSError, ValueError, sqlite3.Error) as error:
                notice = 'Could not create the example report: ' + str(error)
