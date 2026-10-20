"""Shared local input and reference selection for scripted and guided runs."""
import os
from pathlib import Path
from .engine import analyze
from .parser import parse_genome, MAX_BYTES

DEFAULT_REFS = Path(os.environ.get('DNA_RESEARCH_REFS', str(Path.home() / '.cache/dna-research-report/references'))).expanduser()


def read_input(path):
    path = Path(path).expanduser()
    if not path.is_file():
        raise ValueError('Choose an existing, uncompressed DNA text file.')
    if path.stat().st_size > MAX_BYTES:
        raise ValueError('Input exceeds the 64 MiB limit.')
    data = path.read_bytes()
    if len(data) > MAX_BYTES:
        raise ValueError('Input exceeds the 64 MiB limit.')
    return data


def analyze_local(data, build, refs, traits_only=False):
    resolved = parse_genome(data, build).build
    cv = Path(refs) / f'clinvar{resolved}.sqlite'
    gw = Path(refs) / 'gwas.sqlite'
    return analyze(data, build, clinvar=cv if cv.exists() and not traits_only else None,
                   gwas=gw if gw.exists() and not traits_only else None, traits_only=traits_only)
