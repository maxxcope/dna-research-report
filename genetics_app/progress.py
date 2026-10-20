"""TTY progress based on real completed units; line events when redirected."""
import json
import shutil
import sys
import threading
import time
from .report import safe_text


class Progress:
    def __init__(self, stream=None, mode='auto'):
        self.stream = stream or sys.stderr
        self.mode = mode
        self.tty = mode == 'auto' and self.stream.isatty()
        self.started = time.monotonic()
        self.completed = 0
        self.total = 0
        self.name = 'Preparing analysis'
        self.detail = ''
        self._stop = threading.Event()
        self._lock = threading.Lock()
        self._last_plain = None
        self._verbosity = None

    def __enter__(self):
        # HTSlib's informational missing-index/header diagnostics are not progress.
        # Failures still raise exceptions and are retained in stage results.
        try:
            import pysam
            self._verbosity = pysam.set_verbosity(0)
        except ImportError:
            pass
        if self.tty:
            self.stream.write('\x1b[?25l')
            self._thread = threading.Thread(target=self._animate, daemon=True)
            self._thread.start()
        return self

    def event(self, event):
        with self._lock:
            if event['kind'] in {'stage_started', 'stage_finished', 'complete'}:
                self.completed = event['completed']
                self.total = event['total']
                self.name = event.get('name', self.name)
                self.detail = ''
            elif event['kind'] == 'references':
                self.name = 'Verifying reference files'
                self.detail = f"{event['current']}/{event['total']} files"
            else:
                self.detail = safe_text(event.get('message', ''))
        if self.mode == 'json':
            self.stream.write(json.dumps(event) + '\n')
        elif not self.tty:
            if event['kind'] == 'references' and event['current'] not in {1, event['total']}:
                return
            if event['kind'] == 'stage_finished':
                line = f"[{event['completed']}/{event['total']}] {event['name']}: {event['status']}"
            else:
                line = event.get('message') or (self.name + (' - ' + self.detail if self.detail else ''))
            line = safe_text(line)
            if line == self._last_plain:
                return
            self.stream.write(line + '\n')
            self._last_plain = line
        self.stream.flush()

    def _animate(self):
        frame = 0
        while not self._stop.is_set():
            with self._lock:
                count = f'{self.completed}/{self.total} stages' if self.total else 'Preparing'
                width = 18
                filled = int(width * self.completed / self.total) if self.total else 0
                bar = '[' + '#' * filled + '-' * (width - filled) + ']'
                elapsed = int(time.monotonic() - self.started)
                spinner = '|/-' + chr(92)
                line = f"{spinner[frame % 4]} {bar} {count} | {elapsed // 60}:{elapsed % 60:02d}"
                if self.detail:
                    line += ' | ' + self.detail
                line += ' | ' + self.name
            columns = max(20, shutil.get_terminal_size((100, 24)).columns - 1)
            self.stream.write('\r\x1b[2K' + safe_text(line)[:columns])
            self.stream.flush()
            frame += 1
            self._stop.wait(.12)

    def __exit__(self, *exc):
        self._stop.set()
        if self.tty:
            self._thread.join()
            self.stream.write('\r\x1b[2K\x1b[?25h')
            if not exc[0] and self.total:
                self.stream.write(f'[{"#" * 18}] {self.completed}/{self.total} stages processed\n')
            self.stream.flush()
        if self._verbosity is not None:
            import pysam
            pysam.set_verbosity(self._verbosity)
