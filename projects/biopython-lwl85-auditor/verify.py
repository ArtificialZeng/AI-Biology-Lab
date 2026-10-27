"""Read-only packet verification; temporary copied packages receive the repair."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

BASE = Path(__file__).resolve().parent


def compute():
    import Bio
    if Bio.__version__ != '1.88':
        raise RuntimeError('Install the pinned Biopython 1.88 for this release audit')
    source = Path(Bio.__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix='lwl85-audit-') as work:
        work = Path(work)
        shutil.copytree(source, work / 'Bio', ignore=shutil.ignore_patterns('__pycache__'))
        env = os.environ.copy()
        env['PYTHONOPTIMIZE'] = '0'
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        env['PYTHONPATH'] = str(source.parent)
        def run():
            result = subprocess.run([sys.executable, '-B', str(BASE / 'exercise.py')],
                                    env=env, capture_output=True, text=True, check=True)
            return json.loads(result.stdout)
        before = run()
        replacements = [
            ('L[2] * A[1] + L[2] * (A[2] + B[2])', 'L[1] * A[1] + L[2] * (A[2] + B[2])'),
            ('L[2] * B[1] + L[0] * (A[0] + B[0])', 'L[1] * B[1] + L[0] * (A[0] + B[0])'),
        ]
        for relative in ['Align/analysis.py', 'codonalign/codonseq.py']:
            path = work / 'Bio' / relative
            text = path.read_text()
            for old, new in replacements:
                if text.count(old) != 1:
                    raise RuntimeError('Unexpected pinned source: ' + relative)
                text = text.replace(old, new)
            path.write_text(text)
        env['PYTHONPATH'] = str(work)
        after = run()
    for stage, expected in [(before, [False, False, False, True]),
                            (after, [True, True, True, True])]:
        if len(stage) != len(expected):
            raise RuntimeError('Fixture count mismatch')
        for row, agrees in zip(stage, expected):
            if len(row['results']) != 4 or any(r['agrees'] is not agrees for r in row['results']):
                raise RuntimeError('API/order result mismatch: ' + row['name'])
    return {'schema': 'lwl85-release-audit-v1', 'biopython': '1.88',
            'scope': 'Four bounded synthetic standard-code fixtures, both APIs and both orders',
            'float_serialization': '12 decimals; agreement uses 1e-11 relative/1e-12 absolute tolerance',
            'before': before, 'after': after}


def main():
    action = sys.argv[1] if len(sys.argv) == 2 else 'verify'
    if action not in ['verify', 'compute']:
        raise SystemExit('usage: verify.py [verify|compute]')
    if action == 'verify':
        manifest = json.loads((BASE / 'archive_manifest.json').read_text())
        for name, digest in manifest['files'].items():
            if hashlib.sha256((BASE / name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Changed packet input: ' + name)
        actual = compute()
        if actual != json.loads((BASE / 'certificate.json').read_text()):
            raise RuntimeError('Recomputed certificate differs')
        print('PASS: pinned release, temporary repair, 4 fixtures x 2 APIs x 2 orders; packet unchanged')
    else:
        print(json.dumps(compute(), sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
