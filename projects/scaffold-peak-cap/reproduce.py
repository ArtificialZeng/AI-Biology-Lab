#!/usr/bin/env python3
"""Read-only reproduction with assertion-enabled scientific children."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.expanduser().resolve()
    require(work != ROOT and ROOT not in work.parents,
            'Work directory must be outside the project.')
    manifest = json.loads((ROOT / 'manifest.json').read_text())
    for name, expected in manifest['files'].items():
        path = (ROOT / name).resolve()
        require(path.parent == ROOT, 'Manifest contains a nonlocal path.')
        require(path.is_file(), 'Missing manifest file: ' + name)
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        require(actual == expected, 'Manifest mismatch: ' + name)
    work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.pop('PYTHONOPTIMIZE', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    bootstrap = (
        'import sys,runpy; '
        'assert sys.flags.optimize == 0, "Scientific assertions disabled"; '
        'target=sys.argv[1]; sys.argv=sys.argv[1:]; '
        'runpy.run_path(target,run_name="__main__")'
    )
    version = subprocess.run(
        [sys.executable, '-B', '-c',
         'import sys,sympy,json; assert sys.flags.optimize==0; '
         'assert sympy.__version__=="1.13.3"; '
         'print(json.dumps({"python":sys.version.split()[0],'
         '"sympy":sympy.__version__,"optimize":sys.flags.optimize}))'],
        cwd=work, env=env, text=True, capture_output=True)
    (work / 'dependency.stdout').write_text(version.stdout)
    (work / 'dependency.stderr').write_text(version.stderr)
    require(version.returncode == 0, 'Dependency/assertion check failed; see work directory.')
    results = {'dependency': json.loads(version.stdout), 'checks': {}}
    checks = [('coefficient', 'verify_peak_cap.py', []),
              ('sharpness-bisection', 'exact_verifier.py', ['source-baseline.json']),
              ('sharpness-radical', 'exact_verifier.py', ['source-baseline-radical.json']),
              ('one-vacancy', 'check_certificates.py', ['cases.json'])]
    for label, script, inputs in checks:
        command = [sys.executable, '-B', '-c', bootstrap, str(ROOT / script)]
        command += [str(ROOT / item) for item in inputs]
        child = subprocess.run(command, cwd=work, env=env, text=True, capture_output=True)
        (work / (label + '.stdout')).write_text(child.stdout)
        (work / (label + '.stderr')).write_text(child.stderr)
        require(child.returncode == 0, label + ' scientific check failed; see work directory.')
        doc = json.loads(child.stdout)
        if label == 'coefficient':
            require(doc.get('exact_identity_reconstructed') is True and
                    doc.get('integer_coefficients') == 322 and
                    doc.get('all_strictly_negative') is True,
                    'Coefficient reconstruction returned an unexpected result.')
        elif label == 'one-vacancy':
            entries = doc.get('cases', [])
            require(doc.get('status') == 'exact finite certificates valid' and
                    [entry.get('p') for entry in entries] == [2, 3, 4, 5, 6] and
                    all(entry.get('exact_valid') is True and
                        entry.get('certified_disjoint_maximum_intervals') == entry.get('p') and
                        entry.get('certified_signs') == [1, -1] * entry.get('p')
                        for entry in entries),
                    'One-vacancy certificates returned an unexpected result.')
        else:
            require(doc.get('exact_valid') is True and
                    doc.get('signs') == [1, -1, 1, -1] and
                    doc.get('at_least_strict_maxima') == 2,
                    'Sharpness certificate returned an unexpected result.')
        results['checks'][label] = doc
    regenerated = work / 'generated-cases.json'
    command = [sys.executable, '-B', '-c', bootstrap,
               str(ROOT / 'generate_certificates.py'), '--output', str(regenerated)]
    child = subprocess.run(command, cwd=work, env=env, text=True, capture_output=True)
    (work / 'one-vacancy-regeneration.stdout').write_text(child.stdout)
    (work / 'one-vacancy-regeneration.stderr').write_text(child.stderr)
    require(child.returncode == 0, 'One-vacancy regeneration failed.')
    require(regenerated.read_bytes() == (ROOT / 'cases.json').read_bytes(),
            'Regenerated one-vacancy certificates differ from frozen cases.')
    results['checks']['one-vacancy-regeneration'] = {'byte_identical': True}
    results['status'] = 'exact-checks-passed'
    results['scope'] = 'Finite coefficient identity, rational sharpness and one-vacancy case checks; read the universal cap proof and all-p construction proof separately.'
    (work / 'result.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results))

if __name__ == '__main__':
    main()
