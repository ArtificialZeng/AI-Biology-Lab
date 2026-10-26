#!/usr/bin/env python3
"""Build pinned GNU source in temporary trees and check admitted exact-input grids."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ARCHIVE_SHA = '4a1013eebb50f728fc601bdd833b0b2870333c3b3e5a816eeba921d95bec6f15'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def outside(path):
    p = path.resolve()
    require(p != HERE and HERE not in p.parents, 'Keep outputs outside the archived project')
    return p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True, help='official glpk-5.0.tar.gz')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--logs', type=Path)
    args = parser.parse_args()
    require(1 <= args.jobs <= 8, 'Use 1 to 8 build jobs')
    require(sha(args.archive) == ARCHIVE_SHA, 'Unexpected source archive bytes')
    manifest = json.loads((HERE / 'archive_manifest.json').read_text())
    for name, digest in manifest['sha256'].items():
        require(sha(HERE / name) == digest, 'Archived input changed: ' + name)
    output = outside(args.output) if args.output else None
    logs = outside(args.logs) if args.logs else None
    if output:
        require(not output.exists(), 'Refuse to overwrite an existing output receipt')
    if logs:
        logs.mkdir(parents=True, exist_ok=True)
    result = {'source_archive_sha256': ARCHIVE_SHA,
              'patch_sha256': sha(HERE / 'native-guard.patch'),
              'evidence_level': 'finite exact-input solver tests; native algorithm uses floating arithmetic',
              'scope': 'integer singleton lower-rounding branch; existing near-integer tolerance remains',
              'builds': {}}
    with tempfile.TemporaryDirectory(prefix='glpk-presolve-audit-') as folder:
        temp = Path(folder)
        (temp / 'home').mkdir(); (temp / 'tmp').mkdir()
        env = {'PATH': os.defpath, 'HOME': str(temp / 'home'), 'TMPDIR': str(temp / 'tmp'),
               'LANG': 'C', 'LC_ALL': 'C', 'CC': 'cc', 'CFLAGS': '-O2'}
        def run(argv, cwd, label, timeout=180):
            p = subprocess.run([str(x) for x in argv], cwd=cwd, env=env,
                               text=True, capture_output=True, timeout=timeout)
            if logs:
                clean = (p.stdout + p.stderr).replace(str(temp), '<temporary>').replace(str(HERE), '<project>')
                (logs / (label + '.txt')).write_text(clean)
            require(p.returncode == 0, 'Command failed: ' + label + '; request separate logs')
            return p.stdout
        extraction = temp / 'extracted'; extraction.mkdir()
        with tarfile.open(args.archive.resolve(), 'r:gz') as archive:
            archive.extractall(extraction, filter='data')
        source = extraction / 'glpk-5.0'
        require(source.is_dir(), 'Missing pinned source root')
        for label in ('baseline', 'candidate'):
            tree = temp / label
            shutil.copytree(source, tree)
            if label == 'candidate':
                run(['git', 'apply', '--check', HERE / 'native-guard.patch'], tree, 'patch-check')
                run(['git', 'apply', HERE / 'native-guard.patch'], tree, 'patch-apply')
            run(['./configure', '--disable-shared', '--enable-static'], tree, label + '-configure')
            run(['make', '-j' + str(args.jobs)], tree, label + '-build')
            reports = {}
            for harness in ('grid', 'independent_grid'):
                executable = temp / (label + '-' + harness)
                run(['cc', '-O2', '-I' + str(tree / 'src'), HERE / (harness + '.c'),
                     tree / 'src/.libs/libglpk.a', '-lm', '-o', executable], temp,
                    label + '-' + harness + '-compile')
                stdout = run([executable], temp, label + '-' + harness, timeout=120)
                rows = [json.loads(line) for line in stdout.splitlines() if line.startswith('{')]
                if harness == 'grid':
                    require(len(rows) == 1, 'Unexpected primary grid output')
                    report = rows[0]
                    expected = 2600 if label == 'baseline' else 0
                    require(report['cases'] == 72000 and report['expected_feasible'] == 47528,
                            'Primary coverage changed')
                    require(report['wrong_status'] == report['invalid_optimal'] == expected
                            and report['unexpected_return'] == 0, 'Primary grid result changed')
                else:
                    report = next(x for x in rows if x.get('summary'))
                    expected = 12030 if label == 'baseline' else 0
                    require(report['cases'] == 24000 and report['feasible'] == 6134,
                            'Independent coverage changed')
                    require(report['wrong_status'] == report['invalid_optimal'] == expected,
                            'Independent grid result changed')
                    require(report['bad_objective'] == report['bad_recovered_row']
                            == report['unexpected_return'] == 0, 'Feasible recovery regressed')
                    residual = next(x for x in rows if x.get('residual_control') == 1
                                    and x['split'] == 0 and x['presolve'] == 1)
                    require(residual['status'] == 5 and residual['upper_violation'] == 5,
                            'The documented residual changed; reassess scope')
                    reports['known_residual'] = residual
                    controls = [x for x in rows if 'residual_control' in x
                                and (x['split'] == 1 or x['presolve'] == 0)]
                    require(len(controls) == 6 and all(x['status'] == 4 for x in controls),
                            'Presolve-off/split-row controls changed')
                reports[harness] = report
            result['builds'][label] = reports
    rendered = json.dumps(result, indent=2) + '\n'
    if output:
        output.parent.mkdir(parents=True, exist_ok=True); output.write_text(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
