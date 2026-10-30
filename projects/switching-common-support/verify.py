#!/usr/bin/env python3
"""Read-only exact reproduction with external workers and exclusive regeneration."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
MANIFEST = 'archive_manifest.json'
SCHEMA = 'switching-common-support-archive-v1'
CERTS = ('evidence/owner-certificate.json', 'evidence/independent-certificate.json')


class VerificationError(Exception):
    pass


def require(value, message):
    if not value:
        raise VerificationError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def encode(obj):
    return (json.dumps(obj, sort_keys=True, indent=2) + '\n').encode('utf-8')


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON key: ' + key)
        result[key] = value
    return result


def inventory():
    files, directories = [], []
    for item in ROOT.rglob('*'):
        relative = item.relative_to(ROOT).as_posix()
        require(not item.is_symlink(), 'Symlink forbidden: ' + relative)
        if item.is_file():
            if relative != MANIFEST:
                files.append(relative)
        elif item.is_dir():
            directories.append(relative)
        else:
            raise VerificationError('Non-regular archive entry: ' + relative)
    return sorted(files), sorted(directories)


def check_archive():
    require((ROOT / MANIFEST).is_file() and not (ROOT / MANIFEST).is_symlink(),
            'Missing or symlink manifest')
    raw_manifest = (ROOT / MANIFEST).read_bytes()
    document = json.loads(raw_manifest, object_pairs_hook=unique_object)
    require(set(document) == {'schema', 'algorithm', 'files'}, 'Unexpected manifest fields')
    require(document['schema'] == SCHEMA and document['algorithm'] == 'sha256', 'Unknown manifest format')
    expected = document['files']
    require(isinstance(expected, dict) and bool(expected), 'Empty or invalid file map')
    for name, digest in expected.items():
        require(isinstance(name, str) and '\\' not in name and
                PurePosixPath(name).as_posix() == name and
                not PurePosixPath(name).is_absolute() and
                '..' not in PurePosixPath(name).parts and name != MANIFEST,
                'Invalid manifest path')
        require(isinstance(digest, str) and len(digest) == 64 and
                all(c in '0123456789abcdef' for c in digest), 'Invalid SHA256 value')
    files, directories = inventory()
    require(files == sorted(expected), 'Extra or missing archive input')
    required_directories = set()
    for name in expected:
        required_directories.update(parent.as_posix() for parent in PurePosixPath(name).parents
                                    if parent.as_posix() != '.')
    require(directories == sorted(required_directories), 'Extra or missing archive directory')
    data = {name: (ROOT / name).read_bytes() for name in files}
    for name, contents in data.items():
        require(sha(contents) == expected[name], 'Input hash mismatch: ' + name)
    return raw_manifest, data


def compute(data):
    require(sys.version_info >= (3, 12), 'Python 3.12 or newer is required')
    temporary_root = Path(tempfile.gettempdir()).resolve()
    require(temporary_root != ROOT and ROOT not in temporary_root.parents,
            'Temporary worker directory must be outside project')
    computed = {}
    with tempfile.TemporaryDirectory(prefix='switching-common-support-replay-', dir=temporary_root) as folder:
        workspace = Path(folder)
        env = os.environ.copy()
        env.update(PYTHONDONTWRITEBYTECODE='1', PYTHONOPTIMIZE='0')
        env.pop('MATH_JOB_ID', None)
        env.pop('PYTHONPATH', None)
        for label, source, destination in (
            ('owner', 'science/owner_checker.py', CERTS[0]),
            ('fresh', 'science/independent_checker.py', CERTS[1]),
        ):
            job = workspace / label
            job.mkdir()
            script = job / 'check.py'
            script.write_bytes(data[source])
            (job / 'model.json').write_bytes(data['model.json'])
            (job / 'proof.md').write_bytes(data['proof.md'])
            output = job / 'evidence' / 'certificate.json'
            command = [sys.executable, '-B', str(script)]
            process = subprocess.run(command, cwd=job, env=env, capture_output=True, timeout=120)
            require(process.returncode == 0,
                    'Worker failed: ' + label + ': ' + process.stderr.decode('utf-8', errors='replace'))
            if label == 'owner':
                require(output.is_file(), 'Owner worker output is missing')
                result = output.read_bytes()
            else:
                result = process.stdout
            require(isinstance(json.loads(result), dict), 'Worker result is not a JSON object')
            computed[destination] = result
    require(set(computed) == set(CERTS), 'Incomplete computation')
    return computed


def verified_computation():
    manifest, data = check_archive()
    computed = compute(data)
    after_manifest, after_data = check_archive()
    require(after_manifest == manifest and after_data == data, 'Archive changed during verification')
    return manifest, data, computed


def verify():
    manifest, data, computed = verified_computation()
    for name, result in computed.items():
        require(result == data[name], 'Certificate mismatch: ' + name)
    return {'ok': True, 'mode': 'verify', 'archive_modified': False,
            'public_files_checked': len(data), 'worker_python_optimization': 0,
            'certificates': {name: sha(result) for name, result in computed.items()}}


def regenerate(output):
    require(output.is_absolute(), 'Output directory must be absolute')
    require(not output.is_symlink(), 'Output symlink is forbidden')
    target = output.resolve()
    require(target != ROOT and ROOT not in target.parents, 'Output must be outside project')
    require(not target.exists(), 'Output must be a new exclusive directory')
    manifest, data, computed = verified_computation()
    updated = dict(data)
    updated.update(computed)
    new_manifest = encode({'schema': SCHEMA, 'algorithm': 'sha256',
                           'files': {name: sha(contents) for name, contents in updated.items()}})
    target.mkdir(parents=False, exist_ok=False)
    for name, contents in updated.items():
        destination = target / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open('xb') as stream:
            stream.write(contents)
    with (target / MANIFEST).open('xb') as stream:
        stream.write(new_manifest)
    return {'ok': True, 'mode': 'regenerate', 'archive_modified': False,
            'certificates': {name: sha(result) for name, result in computed.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='mode', required=True)
    commands.add_parser('verify')
    regen = commands.add_parser('regenerate')
    regen.add_argument('--output-dir', type=Path, required=True)
    arguments = parser.parse_args()
    try:
        result = verify() if arguments.mode == 'verify' else regenerate(arguments.output_dir)
    except Exception as error:
        print(json.dumps({'ok': False, 'error': str(error)}, sort_keys=True))
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
