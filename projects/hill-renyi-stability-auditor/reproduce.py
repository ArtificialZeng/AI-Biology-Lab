"""Download, verify, build once, audit before/after in fresh subprocesses."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
import venv


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, help='must not already exist; default creates a fresh temporary directory')
    args = parser.parse_args()
    if sys.version_info < (3, 10):
        raise SystemExit('Python >=3.10 required; tested with Python 3.13.5')
    packet = Path(__file__).resolve().parent
    pin = json.loads((packet / 'pin.json').read_text())
    if digest(packet / 'hill-renyi-stable-order.patch') != pin['patch_sha256']:
        raise SystemExit('frozen patch hash mismatch')
    if args.work_dir:
        work = args.work_dir.resolve()
        work.mkdir(parents=True, exist_ok=False)
    else:
        work = Path(tempfile.mkdtemp(prefix='hill-renyi-replay-')).resolve()
    print('Work directory:', work, flush=True)
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    env['PYTHONNOUSERSITE'] = '1'
    commands = []
    def run(command):
        commands.append([str(x) for x in command])
        print('Running:', ' '.join(map(str, command)), flush=True)
        with (work / 'commands.log').open('a') as log:
            log.write('COMMAND ' + repr(commands[-1]) + '\n')
            log.flush()
            subprocess.run(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    archive = work / 'upstream.tar.gz'
    urllib.request.urlretrieve(pin['archive_url'], archive)
    if digest(archive) != pin['archive_sha256']:
        raise SystemExit('upstream archive SHA256 mismatch')
    with tarfile.open(archive) as t:
        # Reject traversal, symlinks and hardlinks; official pinned archive has none.
        for member in t.getmembers():
            target = (work / member.name).resolve()
            if not target.is_relative_to(work) or member.issym() or member.islnk():
                raise SystemExit('unsafe archive member: ' + member.name)
        t.extractall(work)
    source = work / ('scikit-bio-' + pin['commit'])
    def check_hashes(key):
        for name, wanted in pin[key].items():
            if digest(source / name) != wanted:
                raise SystemExit(key + ' mismatch: ' + name)
    check_hashes('before_sha256')
    environment = work / 'venv'
    # Always creates a workdir-owned environment, never installs into the caller.
    venv.EnvBuilder(with_pip=True).create(environment)
    python = environment / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
    uv = shutil.which('uv')
    install = [uv, 'pip', 'install', '--python', str(python)] if uv else [str(python), '-m', 'pip', 'install']
    run(install + ['-r', str(packet / 'requirements.txt')])
    # Complete declared runtime dependencies remain enabled during source install.
    # Frozen build tools are installed above, matching the research build method.
    run(install + ['--no-build-isolation', '-e', str(source)])
    if uv:
        run([uv, 'pip', 'check', '--python', str(python)])
    else:
        run([str(python), '-m', 'pip', 'check'])
    def audit(mode, key):
        run([str(python), '-I', '-B', str(packet / 'audit_api.py'), '--mode', mode, '--source', str(source), '--expected-hash', pin[key]['skbio/diversity/alpha/_base.py'], '--output', str(work / (mode + '.json'))])
    audit('baseline', 'before_sha256')
    # Apply only within the newly downloaded source directory.
    commands.append(['patch', '-p1', '--batch', '-i', str(packet / 'hill-renyi-stable-order.patch')])
    with (work / 'commands.log').open('a') as log:
        subprocess.run(commands[-1], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    check_hashes('after_sha256')
    for cache in (source / 'skbio/diversity/alpha').glob('__pycache__/_base.*.pyc'):
        cache.unlink()
    audit('candidate', 'after_sha256')
    test = source / 'skbio/diversity/alpha/tests/test_base.py'
    run([str(python), '-I', '-B', '-m', 'pytest', '-q', '-o', 'addopts=', str(test)])
    report = {'all_pass': True, 'source_pin': pin['commit'], 'archive_sha256': pin['archive_sha256'], 'patch_sha256': pin['patch_sha256'], 'native_builds': 1, 'method': 'one genuine editable pinned-source native build; baseline public API in fresh process; pure-Python patch; candidate public API in a new process', 'work_dir': str(work), 'commands': commands}
    (work / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS:', work / 'result.json', flush=True)


if __name__ == '__main__':
    main()
