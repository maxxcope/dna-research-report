"""Complete, local report exports with no browser or external assets."""
import json
import gzip
import os
from pathlib import Path
import re
import tempfile
import shutil

FORMATS = {'bundle': None, 'json': 'json', 'markdown': 'md', 'text': 'txt'}
SECTIONS = {
    'knowledge': 'Variant and multi-marker findings', 'conclusions': 'Evidence summaries',
    'clinvar': 'Clinical annotations', 'gwas': 'Study associations', 'traits': 'Curated interpretations',
    'research': 'Additional analyses', 'finding_audit': 'Finding audit', 'qc': 'Input quality',
    'references': 'References', 'pipeline': 'Pipeline completion', 'limits': 'Interpretation limits',
}

def safe_text(value):
    """Do not let source prose or filenames inject terminal control sequences."""
    return re.sub(r'[\x00-\x08\x0b-\x1f\x7f-\x9f]', '', str(value))


def report_lines(report, markdown=True):
    """Walk every field without pagination, truncation, or dropping audit data."""
    yield '# DNA Research Report' if markdown else 'DNA RESEARCH REPORT'
    yield ''
    yield 'Complete retained results, including supporting evidence and unresolved records.'
    yield 'Study associations, source interpretations and statistical estimates retain their original qualifications.'
    yield ''

    def walk(value, level=0):
        prefix = '  ' * level
        if isinstance(value, dict):
            if not value:
                yield prefix + '{}'
            for key, item in value.items():
                label = safe_text(key)
                if isinstance(item, (dict, list)):
                    yield prefix + label + ':'
                    yield from walk(item, level + 1)
                else:
                    lines = safe_text(json.dumps(item, ensure_ascii=False) if not isinstance(item, str) else item).splitlines() or ['']
                    yield prefix + label + ': ' + lines[0]
                    for line in lines[1:]:
                        yield prefix + '  ' + line
        elif isinstance(value, list):
            if not value:
                yield prefix + '[]'
            for number, item in enumerate(value, 1):
                yield prefix + f'[{number}]'
                yield from walk(item, level + 1)
        else:
            yield prefix + safe_text(json.dumps(value, ensure_ascii=False) if not isinstance(value, str) else value)

    for key, value in report.items():
        title = SECTIONS.get(key, key.replace('_', ' ').capitalize())
        yield ('## ' + title) if markdown else title.upper()
        yield 'Field: ' + key
        yield ''
        # Fenced text keeps source prose from becoming Markdown instructions/links.
        if markdown:
            fence = '`' * max(3, 1 + max((len(m[0]) for m in re.finditer(r'`+', json.dumps(value, ensure_ascii=False))), default=0))
            yield fence + 'text'
        yield from walk(value)
        if markdown:
            yield fence
        yield ''


def output_folder(out, create=False):
    """Keep personal results outside the application's distributable folder."""
    out = Path(out).expanduser().resolve()
    project = Path(__file__).resolve().parents[1]
    projects = {project}
    # An installed wheel lives in site-packages, but a user may run it while
    # working in a source checkout. Protect that checkout too, including after
    # renaming it; do not rely on the installation location or folder name.
    for parent in (out, *out.parents):
        metadata = parent / 'pyproject.toml'
        if metadata.is_file():
            text = metadata.read_text(encoding='utf-8')
            section = re.search(r'^\[project\]\s*\n(.*?)(?=^\[|\Z)', text, re.M | re.S)
            if section and re.search(r'''^name\s*=\s*['"]dna-research-report['"]\s*$''', section[1], re.M):
                projects.add(parent)
    if any(out == root or root in out.parents for root in projects):
        raise ValueError('Choose a report folder outside the project to keep personal results out of GitHub.')
    if out.exists() and not out.is_dir():
        raise ValueError('Choose a folder, not a file, for the reports.')
    if create:
        out.mkdir(parents=True, exist_ok=True, mode=0o700)
        with tempfile.TemporaryFile(dir=out):
            pass
    return out


def read_report(path):
    try:
        path = Path(path)
        if path.is_dir():
            path = path / 'full-report.json.gz'
        if path.suffix == '.gz':
            with gzip.open(path, 'rt', encoding='utf-8') as stream:
                report = json.load(stream)
        else:
            report = json.loads(path.read_text(encoding='utf-8'))
    except (gzip.BadGzipFile, EOFError, UnicodeError) as error:
        raise ValueError('Choose a valid JSON report or full-report.json.gz archive.') from error
    except RecursionError as error:
        raise ValueError('The JSON file is nested too deeply to be a supported report.') from error
    if not isinstance(report, dict) or not isinstance(report.get('input_sha256'), str):
        raise ValueError('Choose a complete JSON report generated by this application.')
    if any(not isinstance(report.get(key), list) for key in ('traits', 'clinvar', 'gwas')):
        raise ValueError('Choose a complete JSON report generated by this application.')
    # Saved reports are user-selected files. Reject malformed known sections
    # before exporters or menus try to navigate them; unknown fields stay intact.
    sections = {'': ('traits', 'clinvar', 'gwas', 'exclusions'),
                'knowledge': ('findings', 'unresolved', 'genoset_audit'),
                'finding_audit': ('flagged_findings', 'unresolved_variants', 'genosets'),
                'conclusions': ('cards',)}
    for name, collections in sections.items():
        value = report if not name else report.get(name, {})
        if not isinstance(value, dict):
            raise ValueError('The saved report has an invalid section: ' + name)
        for key in collections:
            if key in value and (not isinstance(value[key], list) or any(not isinstance(item, dict) for item in value[key])):
                raise ValueError('The saved report has invalid records: ' + '.'.join(filter(None, (name, key))))
    if 'research' in report and not isinstance(report['research'], dict):
        raise ValueError('The saved report has an invalid research section.')
    if 'pipeline' in report:
        pipeline = report['pipeline']
        if not isinstance(pipeline, dict) or not isinstance(pipeline.get('status'), str) or not isinstance(pipeline.get('stages'), list) or any(not isinstance(s, dict) or any(not isinstance(s.get(k), str) for k in ('name', 'status')) for s in pipeline['stages']):
            raise ValueError('The saved report has invalid analysis-stage information.')
    if 'knowledge' in report:
        knowledge = report['knowledge']
        if not isinstance(knowledge, dict) or not isinstance(knowledge.get('counts', {}), dict):
            raise ValueError('The saved report has invalid finding counts.')
        counts = knowledge.get('counts')
        if counts and any(not isinstance(counts.get(k), int) or counts[k] < 0 for k in ('variant_summaries', 'genosets', 'unresolved_loci')):
            raise ValueError('The saved report has invalid finding counts.')
    return report


def export_report(report, out, formats=('bundle',)):
    out = output_folder(out)
    formats = tuple(dict.fromkeys(formats))
    if not formats or any(f not in FORMATS for f in formats):
        raise ValueError('Choose bundle, json, markdown or text.')
    targets = {f: out / ('report-bundle' if f == 'bundle' else 'report.' + FORMATS[f]) for f in formats}
    if any(p.exists() or p.is_symlink() for p in targets.values()):
        raise ValueError('An export already exists in this folder. Choose another folder; existing reports are never overwritten.')
    out.mkdir(parents=True, exist_ok=True, mode=0o700)
    with tempfile.TemporaryDirectory(prefix='.export-', dir=out) as tmp:
        staged = {}
        for fmt, target in targets.items():
            path = Path(tmp) / target.name
            if fmt == 'bundle':
                from .bundle import write_bundle
                write_bundle(report, path)
                staged[fmt] = path
                continue
            with path.open('w', encoding='utf-8', newline='\n') as stream:
                if fmt == 'json':
                    json.dump(report, stream, ensure_ascii=True, indent=2)
                    stream.write('\n')
                else:
                    for line in report_lines(report, markdown=fmt == 'markdown'):
                        stream.write(line + '\n')
            os.chmod(path, 0o600)
            staged[fmt] = path
        for fmt, path in staged.items():
            if fmt == 'bundle':
                # Reserve the name atomically, even against an existing empty
                # folder. Only our newly reserved folder is removed on failure.
                targets[fmt].mkdir(mode=0o700)
                try:
                    children = sorted(path.iterdir(), key=lambda p: p.name == 'manifest.json')
                    for child in children:
                        os.rename(child, targets[fmt] / child.name)
                except BaseException:
                    shutil.rmtree(targets[fmt])
                    raise
                continue
            # link is atomic and fails if another process created the target.
            os.link(path, targets[fmt])
    return targets
