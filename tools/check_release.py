"""Check the intended public payload, not ignored local runtime directories.

This catches accidental files/paths; it is not a genomic privacy proof.
"""
from pathlib import Path
import re
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
SKIP={'.git','.local','__pycache__','.venv','build','dist'}
ALLOWED={'.py','.json','.txt','.md','.toml','.yml','.yaml'}
DENIED=['mygenos.push.apply'+'(mygenos,decompressString(']
PRIVATE_PATH=re.compile(r'/(?:Users|home)/[^/\s]+/')
errors=[];count=0
for p in ROOT.rglob('*'):
    rel=p.relative_to(ROOT)
    if any(x in SKIP or x.endswith('.egg-info') for x in rel.parts):continue
    if p.is_symlink():errors.append(f'Symlink: {rel}');continue
    if not p.is_file():continue
    count+=1
    if p.name not in {'.gitignore','LICENSE','start','requirements.lock'} and p.suffix not in ALLOWED:errors.append(f'Unexpected format: {rel}')
    if any(x in {'inputs','analysis','output','references','logs','private','reports','report-bundle','followup'} for x in rel.parts[:-1]):errors.append(f'Non-public directory: {rel}')
    if p.stat().st_size>500_000:errors.append(f'Unexpected large file: {rel}')
    try:text=p.read_text()
    except UnicodeDecodeError:errors.append(f'Binary file: {rel}');continue
    if any(s in text for s in DENIED) or PRIVATE_PATH.search(text):errors.append(f'Private path or original artifact marker: {rel}')
# A .local directory must also stay out of a future Git index.
if (ROOT/'.git').exists():
    tracked=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    for name in tracked:
        if any(x in {'.local','inputs','analysis','output','references','private','reports','report-bundle'} for x in Path(name).parts[:-1]):errors.append(f'Private/runtime path tracked by Git: {name}')
if errors:
    print('\n'.join(errors));sys.exit(1)
print(f'Public payload checks passed: {count} files. Manual publication review is still required.')
