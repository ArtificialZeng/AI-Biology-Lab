"""Shared archive I/O only; the four exact derivations live in separate files."""
from __future__ import annotations

import argparse
import hashlib
import json
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
MANIFEST = 'archive_manifest.json'
BUILDERS = (
    ('evidence/exact-certificate.json', 'verify.py'),
    ('evidence/independent-certificate.json', 'verify_independent.py'),
    ('evidence/four-particle/exact-certificate.json', 'verify_four_particle.py'),
    ('evidence/four-particle/independent-certificate.json', 'verify_four_particle_independent.py'),
)
CERTIFICATES = tuple(certificate for certificate, _ in BUILDERS)


class VerificationError(Exception):
    """Archive or regeneration precondition failed."""


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def encode(value):
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + '\n').encode()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def public_paths():
    paths = []
    for path in ROOT.rglob('*'):
        rel = path.relative_to(ROOT)
        if '__pycache__' in rel.parts or path.suffix in {'.pyc', '.pyo'}:
            continue
        require(not path.is_symlink(), 'Project symlink is not an archive input: ' + str(rel))
        if path.is_file() and rel.as_posix() != MANIFEST:
            paths.append(rel.as_posix())
    return sorted(paths)


def check_archive():
    manifest = json.loads((ROOT / MANIFEST).read_bytes())
    require(manifest.get('schema') == 'ancestral-circulation-archive-v2', 'Unknown manifest schema')
    require(manifest.get('algorithm') == 'sha256', 'Unknown digest algorithm')
    require(sorted(manifest.get('files', {})) == public_paths(), 'Manifest must cover every public project input')
    for name, expected in manifest['files'].items():
        require(sha256((ROOT / name).read_bytes()) == expected, 'Input hash mismatch: ' + name)
    return manifest


def verify(builder, certificate, implementation):
    manifest = check_archive()
    recomputed = encode(builder())
    archived = (ROOT / certificate).read_bytes()
    require(recomputed == archived, 'Recomputed certificate differs from archive: ' + certificate)
    return {'ok': True, 'mode': 'verify', 'implementation': implementation,
            'archive_modified': False, 'public_files_checked': len(manifest['files']),
            'certificate_sha256': sha256(archived)}


def load_builder(relative):
    """Load an entry point by its local path; mathematical implementations stay separate."""
    name = '_ancestral_builder_' + relative.replace('/', '_').replace('.', '_')
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    require(spec is not None and spec.loader is not None, 'Cannot load builder: ' + relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_certificate


def regenerate(directory):
    require(directory.is_absolute(), 'Regeneration output must be an absolute path')
    output = directory.expanduser().resolve()
    require(output != ROOT and ROOT not in output.parents, 'Regeneration output must be outside the project')
    require(not output.exists(), 'Regeneration output must be a new exclusive directory')
    # Explicit regeneration independently recomputes all four certificates.
    candidates = {certificate: encode(load_builder(relative)()) for certificate, relative in BUILDERS}
    paths = public_paths()
    require(all(name in paths for name in CERTIFICATES), 'All four archived certificates are required')
    manifest = {'schema': 'ancestral-circulation-archive-v2', 'algorithm': 'sha256',
                'files': {name: sha256(candidates.get(name, (ROOT / name).read_bytes())) for name in paths}}
    output.mkdir(parents=False, exist_ok=False)
    for name, data in candidates.items():
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    (output / MANIFEST).write_bytes(encode(manifest))
    return {'ok': True, 'mode': 'regenerate', 'archive_modified': False,
            'exclusive_output': str(output),
            'certificate_sha256': {name: sha256(data) for name, data in candidates.items()}}


def run(builder, certificate, implementation):
    parser = argparse.ArgumentParser(description='Read-only exact verification; explicit external regeneration only.')
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('verify')
    sub.add_parser('compute')
    regen = sub.add_parser('regenerate')
    regen.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == 'compute':
            print(encode(builder()).decode(), end='')
            return 0
        result = verify(builder, certificate, implementation) if args.command == 'verify' else regenerate(args.output_dir)
    except Exception as exc:
        print(json.dumps({'ok': False, 'error': str(exc)}, sort_keys=True))
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0
