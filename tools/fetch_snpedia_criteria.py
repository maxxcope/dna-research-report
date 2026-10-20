"""Download public genoset criteria for the research audit, without DNA input.

Provide a local SNP Browser schema reference snapshot. This reads every genoset
ID from the public database, not a person's marker list. SNPedia content remains
subject to CC BY-NC-SA 3.0 US. Keep this cache outside the release directory.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import re
from pathlib import Path
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request

API = 'https://bots.snpedia.com/api.php'


def download(database, output):
    with sqlite3.connect(database.resolve().as_uri() + '?mode=ro', uri=True) as connection:
        titles = [row[0].capitalize() + '/criteria' for row in
                  connection.execute('SELECT id FROM genosets ORDER BY id')]
    if not titles or len(titles) > 10_000:
        raise ValueError('Unexpected public genoset inventory size')
    if any(not re.fullmatch(r'Gs\d+/criteria', title) for title in titles):
        raise ValueError('Unexpected genoset identifier')
    source_hash = hashlib.sha256(database.read_bytes()).hexdigest()
    pages = []
    if output.exists():
        previous = json.loads(output.read_text())
        if previous.get('database_sha256') != source_hash:
            raise ValueError('Existing checkpoint belongs to another source snapshot')
        pages = previous['pages']
    completed = {page['title'] for page in pages}
    pending = [title for title in titles if title not in completed]
    for start in range(0, len(pending), 40):
        batch = pending[start:start + 40]
        parameters = dict(action='query', prop='revisions', rvprop='ids|timestamp|sha1|content',
                          rvslots='main', titles='|'.join(batch), format='json', formatversion='2')
        request = urllib.request.Request(API + '?' + urllib.parse.urlencode(parameters),
            headers={'User-Agent': 'GeneticsParityResearch/0.1 (public genoset rule audit)'})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request, timeout=45) as response:
                    data = json.load(response)
                break
            except urllib.error.HTTPError as error:
                if error.code not in {429, 502, 503, 504} or attempt == 2:
                    raise
                retry = error.headers.get('Retry-After', '')
                delay = float(retry) if retry.isdigit() else 0
                time.sleep(max(delay, 10 * (attempt + 1)))
        if 'error' in data or 'warnings' in data:
            raise ValueError('SNPedia returned an API error or warning')
        fetched = data['query']['pages']
        if {page['title'] for page in fetched} != set(batch):
            raise ValueError('Response did not contain the requested public rule titles')
        for page in fetched:
            if 'revisions' not in page:
                if 'missing' not in page:
                    raise ValueError('Rule neither readable nor explicitly missing')
                continue
            revision = page['revisions'][0]
            content = revision['slots']['main']['content']
            if hashlib.sha1(content.encode()).hexdigest() != revision['sha1']:
                raise ValueError('Public rule SHA1 mismatch')
        pages.extend(fetched)
        record = dict(source=API, license='CC BY-NC-SA 3.0 US',
                      license_url='https://creativecommons.org/licenses/by-nc-sa/3.0/us/',
                      database_sha256=source_hash, retrieved_at=datetime.now(timezone.utc).isoformat(),
                      complete=len(pages) == len(titles), pages=pages)
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_suffix(output.suffix + '.tmp')
        temporary.write_text(json.dumps(record, indent=2) + '\n')
        temporary.replace(output)
        print(f'Fetched {len(pages)} of {len(titles)} public rule pages', flush=True)
        if start + 40 < len(pending):
            time.sleep(3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.resolve().is_relative_to(Path(__file__).resolve().parents[1]):
        parser.error('Keep the third-party data cache outside the release folder')
    if args.out.resolve() == args.database.resolve():
        parser.error('Output must not overwrite the reference database')
    download(args.database, args.out)


if __name__ == '__main__':
    main()
