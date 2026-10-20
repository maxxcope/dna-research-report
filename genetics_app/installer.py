"""First-run installation. Only public references are downloaded here."""
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from .downloads import atomic_json, download, sha256

ROOT = Path(__file__).resolve().parents[1]
MIN_FREE = 10 * 1024 ** 3


def platform_key():
    machine = {'aarch64': 'arm64', 'AMD64': 'x86_64'}.get(platform.machine(), platform.machine())
    key = platform.system().lower() + '-' + machine
    if key not in {'darwin-arm64', 'darwin-x86_64', 'linux-x86_64'}:
        raise ValueError('Automatic setup supports macOS Apple Silicon/Intel and Linux x86_64 with glibc. See docs/REFERENCE_SETUP.md.')
    if key.startswith('linux') and platform.libc_ver()[0] != 'glibc':
        raise ValueError('Automatic Linux setup requires glibc; Alpine/musl is not supported.')
    return key


def external_folder(value, project=ROOT):
    path = Path(value).expanduser().resolve()
    if path == project or project in path.parents:
        raise ValueError('Keep application data and references outside this project. Set DNA_RESEARCH_HOME or DNA_RESEARCH_REFS to another folder.')
    return path


@contextmanager
def setup_lock(home):
    # flock is released by the OS even after SIGKILL. Never delete the lock inode:
    # another waiting launcher may already have it open.
    with (home / 'setup.lock').open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Setup is already running in another terminal. Let it finish, then run ./start again.') from None
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def run(command, log, cwd=ROOT):
    from .report import safe_text
    process = subprocess.Popen([str(x) for x in command], cwd=cwd, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, start_new_session=True)
    try:
        with log.open('a', encoding='utf-8') as output:
            for line in iter(process.stdout.readline, b''):
                text = safe_text(line.decode('utf-8', errors='replace')).rstrip()
                print(text, flush=True)
                output.write(text + '\n')
                output.flush()
        if process.wait():
            raise ValueError('A setup step failed. Details: ' + str(log))
    finally:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
        process.stdout.close()


def locked_versions(lock):
    return dict(re.findall(r'^([a-zA-Z0-9_-]+)==([^\s\\]+)', lock.read_text(), re.M))


def environment_ready(python, versions):
    if not python.is_file():
        return False
    check = """import importlib.metadata as m,json,sys
expected=json.loads(sys.argv[1])
assert all(m.version(k)==v for k,v in expected.items())
import numpy,sklearn,pysam,scipy
"""
    try:
        return subprocess.run([str(python), '-I', '-c', check, json.dumps(versions)],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=60).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def ensure_environment(home, uv, log):
    lock = ROOT / 'requirements.lock'
    versions = locked_versions(lock)
    if not versions:
        raise ValueError('The dependency lock is missing or empty.')
    python = home / 'environment/bin/python'
    if not environment_ready(python, versions):
        print('Installing the analysis environment. Packages are checked against the dependency lock.', flush=True)
        if not python.exists() or not environment_ready(python, {}):
            # A moved or interrupted venv is disposable. Do not touch references.
            if python.parent.parent.exists():
                shutil.rmtree(python.parent.parent)
            run([uv, 'venv', '--python', sys.executable, home / 'environment'], log)
        run([uv, 'pip', 'sync', '--python', python, '--require-hashes', '--only-binary', ':all:', lock], log)
        if not environment_ready(python, versions):
            raise ValueError('Installed packages did not pass their import/version checks. See ' + str(log))
    return python


def executable_works(path, name):
    try:
        result = subprocess.run([str(path), '--version' if name == 'plink' else '-version'],
                                capture_output=True, text=True, timeout=30)
        text = result.stdout + result.stderr
        if result.returncode:
            return False
        if name == 'plink':
            return bool(re.search(r'PLINK v1\.(?:9\d*|90)', text))
        match = re.search(r'version "(?:1\.)?(\d+)', text)
        return bool(match and int(match[1]) >= 8)
    except (OSError, subprocess.TimeoutExpired):
        return False


def extract_tool(archive, destination, name):
    if name == 'plink':
        with zipfile.ZipFile(archive) as source:
            member = source.getinfo('plink')
            if member.is_dir() or member.file_size > 100 * 1024 ** 2:
                raise ValueError('Unexpected PLINK archive.')
            with source.open(member) as src, (destination / 'plink').open('wb') as target:
                shutil.copyfileobj(src, target)
        path = destination / 'plink'
        path.chmod(0o700)
    else:
        # Python 3.12's data filter rejects traversal and escaping links, while
        # allowing the internal library symlinks needed by a JRE.
        with tarfile.open(archive) as source:
            source.extractall(destination, filter='data')
        found = list(destination.rglob('bin/java'))
        if len(found) != 1:
            raise ValueError('The Java archive did not contain one runtime.')
        path = found[0]
    return path


def ensure_tool(home, name, spec):
    override = os.environ.get('DNA_RESEARCH_' + name.upper())
    if override:
        path = Path(shutil.which(override) or override).expanduser().resolve()
        if not executable_works(path, name):
            raise ValueError(f'The configured {name} executable could not run: {path}')
        return path
    folder = home / 'tools' / (name + '-' + spec['sha256'][:16])
    folder.parent.mkdir(parents=True, exist_ok=True)
    for stale in folder.parent.glob('.install-' + name + '-*'):
        if stale.is_dir() and not stale.is_symlink():
            shutil.rmtree(stale)
    record = folder / 'installed.json'
    if record.exists():
        try:
            state = json.loads(record.read_text())
            path = folder / state['executable']
            intact = bool(state['files']) and state['archive_sha256'] == spec['sha256'] and all(
                (folder / r['file']).is_file() and sha256(folder / r['file']) == r['sha256'] for r in state['files'])
            if intact and executable_works(path, name):
                return path
        except (ValueError, KeyError, TypeError, AttributeError, OSError):
            pass
        raise ValueError(f'The installed {name} files changed. Move {folder} aside and run ./start again.')
    archive = download(spec['url'], home / 'downloads' / Path(spec['url']).name,
                       expected=spec['sha256'], size=spec['bytes'])
    folder.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.install-' + name + '-', dir=folder.parent) as temp:
        temporary = Path(temp)
        path = extract_tool(archive, temporary, name)
        if not executable_works(path, name):
            raise ValueError(f'The downloaded {name} cannot run on this system. See docs/REFERENCE_SETUP.md for supported systems and overrides.')
        state = dict(executable=str(path.relative_to(temporary)), archive_sha256=spec['sha256'],
                     files=[dict(file=str(p.relative_to(temporary)), sha256=sha256(p))
                            for p in sorted(temporary.rglob('*')) if p.is_file()])
        atomic_json(temporary / 'installed.json', state)
        # A previous interrupted install never has a completed record.
        if folder.exists():
            shutil.rmtree(folder)
        temporary.rename(folder)
    return folder / state['executable']


def profile_files(config):
    """Reject partial/malformed profiles, including untracked configured files."""
    sections = dict(quality={'clusters', 'badalleles'}, frequencies={'prefix', 'labels'},
                    paternal={'tree', 'chain'}, maternal={'directory', 'java', 'tree'},
                    imputation={'windows', 'java', 'beagle', 'score'}, phase={'window', 'java', 'beagle', 'sites'})
    required = {'knowledge', 'clinvar', 'gwas', 'ancestry', 'scores', 'reference_files', *sections}
    if not isinstance(config, dict) or not all(config.get(k) for k in required):
        raise ValueError('The reference profile is incomplete. Choose a new cache folder; setup will not silently reduce the report.')
    for section, fields in sections.items():
        if not isinstance(config[section], dict) or not fields <= config[section].keys():
            raise ValueError('Incomplete reference profile section: ' + section)
    for section, fields in [('ancestry', {'prefix', 'labels', 'plink'}), ('scores', {'id', 'path', 'build'})]:
        if not isinstance(config[section], list) or any(not isinstance(r, dict) or not fields <= r.keys() for r in config[section]):
            raise ValueError('Malformed reference profile section: ' + section)
    records = config['reference_files']
    if not isinstance(records, list) or any(not isinstance(r, dict) or not isinstance(r.get('file'), str)
           or not re.fullmatch(r'[0-9a-f]{64}', str(r.get('sha256', ''))) for r in records):
        raise ValueError('Malformed reference file inventory. Choose a new cache folder.')
    tracked = {Path(r['file']).resolve() for r in records}
    path_keys = {'clusters', 'badalleles', 'labels', 'plink', 'tree', 'chain', 'java', 'beagle', 'score', 'file', 'map', 'path', 'knowledge', 'clinvar', 'gwas'}
    def check(value, key=''):
        if isinstance(value, dict):
            for name, item in value.items():
                if name != 'reference_files':
                    check(item, name)
        elif isinstance(value, list):
            for item in value:
                check(item, key)
        elif key == 'prefix':
            for suffix in ['bed', 'bim', 'fam']:
                check(str(value) + '.' + suffix, 'file')
        elif key in path_keys:
            if key == 'tree' and isinstance(value, str) and '@' in value:
                return
            if not isinstance(value, str) or Path(value).resolve() not in tracked:
                raise ValueError('Configured reference is missing from the hash inventory: ' + key)
    check(config)
    return records


def validate_profile(refs, marker):
    from .pipeline import load_config
    from .references import metadata
    profile = refs / 'pipeline.json'
    if not profile.exists():
        return False
    try:
        config = load_config(profile)
        files = profile_files(config)
    except (TypeError, KeyError, AttributeError) as error:
        raise ValueError('Malformed reference profile. Choose a new cache folder.') from error
    if any(not Path(r['file']).is_file() for r in files):
        return False
    stamps = [[r['file'], Path(r['file']).stat().st_size, Path(r['file']).stat().st_mtime_ns] for r in files]
    identity = dict(profile_sha256=sha256(profile), files=stamps)
    previous = None
    if marker.exists():
        try:
            previous = json.loads(marker.read_text())
        except (ValueError, UnicodeError):
            pass
    if previous != identity:
        print('Checking prepared reference files. This may take a few minutes.', flush=True)
        for index, record in enumerate(files, 1):
            if sha256(record['file']) != record['sha256']:
                raise ValueError('Reference checksum changed: ' + Path(record['file']).name + '. Move the damaged file aside or choose a new DNA_RESEARCH_REFS folder.')
            if index == 1 or index % 25 == 0 or index == len(files):
                print(f'Reference files checked: {index}/{len(files)}', flush=True)
        if metadata(config['clinvar']).get('build') != '37':
            raise ValueError('Automatic full setup expects a GRCh37 reference profile.')
    tools = [(p['plink'], 'plink') for p in config['ancestry']]
    tools.extend((config[section]['java'], 'java') for section in ['maternal', 'imputation', 'phase'])
    if not all(executable_works(p, name) for p, name in tools):
        # Re-prepare to update paths to working managed tools.
        return False
    if previous != identity:
        atomic_json(marker, identity)
    return True


def setup(home, refs, uv, log):
    python = ensure_environment(home, uv, log)
    marker = home / 'references-ready.json'
    if not validate_profile(refs, marker):
        if shutil.disk_usage(refs if refs.exists() else home).free < MIN_FREE:
            raise ValueError('Full reference preparation needs at least 10 GiB of free disk space. Free space and run ./start again.')
        print('Preparing full public references (GRCh37). Allow tens of minutes on the first run.\nCompleted files are reused after interruption. No personal DNA is read or uploaded.', flush=True)
        specs = json.loads((ROOT / 'genetics_app/data/install_tools.json').read_text())[platform_key()]
        plink = ensure_tool(home, 'plink', specs['plink'])
        java = ensure_tool(home, 'java', specs['java'])
        run([python, '-u', '-m', 'genetics_app.setup', '--refs', refs, '--plink', plink, '--java', java], log)
        if not validate_profile(refs, marker):
            raise ValueError('Reference preparation finished without a valid full profile. See ' + str(log))
    return python


def launch(python, args, refs, limited=False):
    environment = os.environ.copy()
    environment.pop('PYTHONPATH', None)
    environment.pop('PYTHONHOME', None)
    environment['DNA_RESEARCH_REFS'] = str(refs)
    environment['PYTHONUNBUFFERED'] = '1'
    os.chdir(ROOT)
    os.execve(str(python), [str(python), '-I', str(ROOT / 'launch.py'), '--internal-app', *(['--exports-only'] if limited else []), *args], environment)


def _main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if args[:1] == ['--internal-app']:
        from .__main__ import main as application
        return application(args[1:])
    try:
        platform_key()
        home = external_folder(os.environ.get('DNA_RESEARCH_HOME', Path.home() / '.local/share/dna-research-report'))
        refs = external_folder(os.environ.get('DNA_RESEARCH_REFS', home / 'references'))
        home.mkdir(parents=True, exist_ok=True, mode=0o700)
        uv = home / ('uv-0.13.4-' + {'darwin-arm64': 'aarch64-apple-darwin', 'darwin-x86_64': 'x86_64-apple-darwin', 'linux-x86_64': 'x86_64-unknown-linux-gnu'}[platform_key()])
        log = home / 'setup.log'
        app_only = args[:1] == ['--app-only']
        setup_only = args[:1] == ['--setup-only']
        if app_only or setup_only:
            args.pop(0)
        while True:
            try:
                if app_only:
                    # Exports and the small example use only the standard library.
                    # They remain available even if package installation failed.
                    python = Path(sys.executable)
                    break
                with setup_lock(home):
                    python = setup(home, refs, uv, log)
                break
            except (OSError, ValueError, subprocess.SubprocessError) as error:
                from .terminal import Terminal, QuitFlow
                terminal = Terminal()
                terminal.say('\nSetup could not finish: ' + str(error))
                terminal.say('Completed downloads are kept. Run ./start again to retry.\nSetup log: ' + str(log))
                if not terminal.arrows or setup_only:
                    return 1
                try:
                    choice = terminal.choose('What would you like to do?', [('retry', 'Retry setup'), ('app', 'Open the app for exports or an example')], back=False)
                except QuitFlow:
                    return 0
                app_only = choice == 'app'
        if setup_only:
            print('Setup complete. Run ./start to open DNA Research Report.', flush=True)
            return 0
        launch(python, args, refs, limited=app_only)
    except KeyboardInterrupt:
        print('\nSetup cancelled. Completed downloads are kept; run ./start to continue.', file=sys.stderr)
        return 130
    except (OSError, ValueError) as error:
        print('Setup could not continue: ' + str(error), file=sys.stderr)
        return 1


def main(argv=None):
    def cancelled(signum, frame):
        raise KeyboardInterrupt()
    # Closing a terminal must stop a reference-preparation child, too.
    previous = {signum: signal.signal(signum, cancelled) for signum in (signal.SIGTERM, signal.SIGHUP)}
    try:
        return _main(argv)
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)
