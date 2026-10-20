"""Atomic public downloads with bounded retries and validator-guarded resume."""
import hashlib
import http.client
import json
from pathlib import Path
import re
import sys
import time
import urllib.error
import urllib.request


class RestartDownload(ValueError):
    """A server/proxy did not honor the guarded byte range."""


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def byte_progress(name, current, total):
    def size(n):
        if n < 1024:return f'{n} B'
        if n < 1048576:return f'{n / 1024:.1f} KiB'
        return f'{n / 1048576:.1f} MiB'
    suffix = ' / ' + size(total) if total else ''
    bar = ''
    if total:
        completed = min(24, current * 24 // total)
        bar = '[' + '=' * completed + ' ' * (24 - completed) + '] '
    print(f'{bar}{name}: {size(current)}{suffix}', flush=True)


def download(url, path, expected=None, size=None, progress=byte_progress, attempts=3):
    """Never publish partial files. Resume only with a stable server validator."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if not path.is_file() or (expected and sha256(path) != expected):
            raise ValueError('Cached file checksum changed: ' + path.name + '. Move this file aside before retrying setup.')
        if size is not None and path.stat().st_size != size:
            raise ValueError('Cached file has an unexpected size: ' + path.name)
        return path
    part = path.with_name(path.name + '.part')
    checkpoint = path.with_name(path.name + '.part.json')
    for attempt in range(attempts):
        state = {}
        if checkpoint.exists():
            try:
                value = json.loads(checkpoint.read_text())
                if isinstance(value, dict) and isinstance(value.get('validator'), str):
                    state = value
            except (ValueError, UnicodeError):
                pass
        offset = part.stat().st_size if part.exists() else 0
        if state.get('url') != url or not state.get('validator'):
            offset = 0
        headers = {'User-Agent': 'DNAResearchReport/1.0 public-reference-setup', 'Accept-Encoding': 'identity'}
        if offset:
            headers.update({'Range': f'bytes={offset}-', 'If-Range': state['validator']})
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=60) as response:
                status = response.status
                validator = response.headers.get('ETag', '')
                if validator.startswith('W/'):
                    validator = ''
                validator = validator or response.headers.get('Last-Modified', '')
                length = response.headers.get('Content-Length')
                length = int(length) if length and length.isdigit() else None
                if status == 206:
                    match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
                    if not match or int(match[1]) != offset or (offset and validator != state.get('validator')):
                        raise RestartDownload('Server returned an incompatible partial download: ' + path.name)
                    total = int(match[3])
                    if int(match[2]) != total - 1:
                        raise RestartDownload('Server returned an incomplete byte range: ' + path.name)
                elif status == 200:
                    offset = 0
                    total = length
                else:
                    raise ValueError('Unexpected download response: ' + str(status))
                if size is not None and total is not None and total != size:
                    raise ValueError('Download size differs from the pinned release: ' + path.name)
                atomic_json(checkpoint, dict(url=url, validator=validator, total=total))
                current = offset
                last = 0.0
                with part.open('ab' if offset else 'wb') as target:
                    while chunk := response.read(1024 * 1024):
                        target.write(chunk)
                        current += len(chunk)
                        if size is not None and current > size:
                            raise ValueError('Download exceeds the pinned release size: ' + path.name)
                        if time.monotonic() - last >= 1:
                            progress(path.name, current, total)
                            last = time.monotonic()
                if total is not None and current != total:
                    raise http.client.IncompleteRead(b'', total - current)
                if size is not None and current != size:
                    raise ValueError('Incomplete pinned release: ' + path.name)
                if not current:
                    raise ValueError('Empty download: ' + path.name)
                if expected and sha256(part) != expected:
                    part.unlink(missing_ok=True)
                    checkpoint.unlink(missing_ok=True)
                    raise ValueError('Downloaded file checksum did not match: ' + path.name)
                part.replace(path)
                checkpoint.unlink(missing_ok=True)
                progress(path.name, current, total)
                return path
        except RestartDownload:
            part.unlink(missing_ok=True)
            checkpoint.unlink(missing_ok=True)
            if attempt == attempts - 1:
                raise
        except urllib.error.HTTPError as error:
            if error.code == 416 and offset:
                # An interrupted write may already contain the final byte. A fresh
                # request is safer than accepting an unverified partial reference.
                part.unlink(missing_ok=True)
                checkpoint.unlink(missing_ok=True)
            elif error.code not in {408, 429, 500, 502, 503, 504}:
                raise
            if attempt == attempts - 1:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException):
            if attempt == attempts - 1:
                raise
        print(f'Download interrupted; retrying {path.name} ({attempt + 2}/{attempts}).', file=sys.stderr, flush=True)
        time.sleep(min(2 ** attempt, 4))
    raise ValueError('Download did not complete: ' + path.name)
