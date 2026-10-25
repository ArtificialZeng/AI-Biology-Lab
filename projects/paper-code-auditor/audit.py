#!/usr/bin/env python3
"""Reproduce pinned software-entrypoint defects and verify two local fixes.

Runs selected data-handling AST nodes and reconstructed loop control only.
Never imports model libraries or launches upstream training commands.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
UPSTREAM = ROOT / 'fixtures/upstream'
WRAPPER = 'scripts/run_ESM2ne_cellular.sh'
TRAINER = 'scripts/LitESM2ne_Lora_trainer.py'
ASSAY = ROOT / 'fixtures/data/IF1_ECOLI_identifiers.csv'


def integrity():
    manifest = json.loads((ROOT/'evidence/input_manifest.json').read_text())
    for name, expected in manifest['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Input hash mismatch: {name}')
    return len(manifest['sha256'])


def functions(text, names):
    nodes = [node for node in ast.parse(text).body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in nodes} == set(names)
    namespace = {'np': np, 'pd': pd, 'hashlib': hashlib}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'selected_pinned_functions', 'exec'), namespace)
    return namespace


def setup_function(text, seed):
    tree = ast.parse(text)
    klass = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'MyDataModule')
    setup = next(n for n in klass.body if isinstance(n, ast.FunctionDef) and n.name == 'setup')
    scope = functions(text, {'split_data'})
    scope.update(args=SimpleNamespace(data=str(ASSAY), seed=seed), MyDataset=lambda frame: frame)
    exec(compile(ast.Module(body=[setup], type_ignores=[]), 'selected_pinned_setup', 'exec'), scope)
    return scope['setup']


def dispatch_count(text):
    """Reconstruct reviewed loop heads, replacing the entire body by a counter."""
    import re
    datasets = shlex.split(text.split('cellular=(', 1)[1].split(')', 1)[0], comments=True)
    strategies = shlex.split(text.split('strategies=(', 1)[1].split(')', 1)[0])
    variable = re.search(r'for dts in "\$\{([a-zA-Z_][a-zA-Z_0-9]*)\[@\]\}"', text).group(1)
    replicates = re.search(r'for rep in ([0-9 ]+)\ndo', text).group(1).split()
    assert len(datasets) == 33 and strategies == ['pool_split', 'site_split']
    assert replicates == ['1', '2', '3']
    assert variable in {'viral', 'cellular'} and 'viral=(' not in text
    control = ('python() { printf "CALL\\n"; }\n'
               + 'cellular=(' + ' '.join(shlex.quote(n) for n in datasets) + ')\n'
               + 'strategies=(' + ' '.join(shlex.quote(n) for n in strategies) + ')\n'
               + 'for rep in 1 2 3; do\n'
               + ' for strategy in "${strategies[@]}"; do\n'
               + '  for dts in "${' + variable + '[@]}"; do python; done\n'
               + ' done\ndone\n')
    bash = shutil.which('bash')
    if not bash:
        raise RuntimeError('Bash is required for the isolated control-flow check')
    result = subprocess.run([bash, '--noprofile', '--norc', '-c', control],
                            env={'PATH': '/usr/bin:/bin'}, capture_output=True, text=True, check=True)
    return {'exit_code': result.returncode, 'counter_calls': result.stdout.count('CALL\n'),
            'datasets': len(datasets), 'strategies': len(strategies), 'replicates': len(replicates)}


def check():
    file_count = integrity()
    original_text = (UPSTREAM/TRAINER).read_text()
    lasso = functions((UPSTREAM/'scripts/lassoCV_SS.py').read_text(), {'make_seed', 'split_data'})
    original_helpers = functions(original_text, {'split_data'})
    data = pd.read_csv(ASSAY).query("mutant != 'WT'").copy()
    assert len(data) == 1439
    inventory = json.loads((ROOT/'evidence/upstream_inventory.json').read_text())
    missing = 'scripts/LitESM2ne_partial_trainer.py'
    assert inventory['truncated'] is False and missing not in inventory['files']
    assert missing in (UPSTREAM/WRAPPER).read_text()
    assert missing in (UPSTREAM/'scripts/run_ESM2ne_viral.sh').read_text()
    original_dispatch = dispatch_count((UPSTREAM/WRAPPER).read_text())
    assert original_dispatch['counter_calls'] == 0
    seed = lasso['make_seed']('IF1_ECOLI', 1)
    try:
        setup_function(original_text, seed)(SimpleNamespace(split_strategy='site_split'))
    except TypeError as error:
        message = str(error)
        assert 'train_pct' in message
    else:
        raise AssertionError('Expected the pinned unsupported-keyword failure')

    with tempfile.TemporaryDirectory(prefix='paper-code-check-') as temporary:
        work = Path(temporary)
        shutil.copytree(UPSTREAM, work, dirs_exist_ok=True)
        patch = ROOT/'patches/local-fixes.patch'
        subprocess.run(['git', 'apply', '--check', str(patch)], cwd=work, check=True, capture_output=True)
        subprocess.run(['git', 'apply', str(patch)], cwd=work, check=True, capture_output=True)
        patched_text = (work/TRAINER).read_text()
        fixed_dispatch = dispatch_count((work/WRAPPER).read_text())
        assert fixed_dispatch['counter_calls'] == 198
        # The partial trainer is deliberately not replaced by a different method.
        assert missing in (work/WRAPPER).read_text()
        rows = []
        for replicate in [1, 2, 3]:
            seed = lasso['make_seed']('IF1_ECOLI', replicate)
            a, b = lasso['split_data'](data, seed)
            c, d = original_helpers['split_data'](data, seed)
            assert set(b.ID) == set(d.ID)
            obj = SimpleNamespace(split_strategy='site_split')
            setup_function(patched_text, seed)(obj)
            train, valid = obj.train_dataset, obj.val_dataset
            assert set(valid.ID) == set(b.ID)
            assert len(train) + len(valid) == len(data)
            assert set(train.site).isdisjoint(valid.site)
            assert set(train.ID).isdisjoint(valid.ID)
            rows.append({'replicate': replicate, 'seed': seed, 'train_rows': len(train),
                         'validation_rows': len(valid), 'sites_disjoint': True,
                         'all_rows_preserved': True, 'matches_Lasso_fold': True,
                         'original_Lasso_LoRA_helpers_agree': True})
        before = ast.parse(original_text)
        after = ast.parse(patched_text)
        original_split = next(n for n in before.body if isinstance(n, ast.FunctionDef) and n.name == 'split_data')
        patched_split = next(n for n in after.body if isinstance(n, ast.FunctionDef) and n.name == 'split_data')
        assert ast.dump(original_split) == ast.dump(patched_split)
    restoration = json.loads((ROOT/'evidence/restoration_provenance.json').read_text())
    with tempfile.TemporaryDirectory(prefix='paper-code-restore-') as temporary:
        work = Path(temporary)
        shutil.copytree(UPSTREAM, work, dirs_exist_ok=True)
        patch = ROOT/'patches/reproduction-fixes.patch'
        subprocess.run(['git', 'apply', '--check', str(patch)], cwd=work, check=True, capture_output=True)
        subprocess.run(['git', 'apply', str(patch)], cwd=work, check=True, capture_output=True)
        for name, metadata in restoration['source_files'].items():
            historical = (ROOT/'fixtures/historical'/name).read_bytes()
            assert hashlib.sha256(historical).hexdigest() == metadata['sha256']
            expected = historical
            if name.endswith('LitESM2ne_partial_trainer.py'):
                expected = historical.replace(b'devices=[1],', b'devices=[0],')
            assert (work/name).read_bytes() == expected
        historical_text = (work/missing).read_text()
        assert 'from scripts.ESM2ne_partial import LoadFromPretrained' in historical_text
        model_tree = ast.parse((work/'scripts/ESM2ne_partial.py').read_text())
        assert any(isinstance(n, ast.ClassDef) and n.name == 'LoadFromPretrained' for n in model_tree.body)
        for wrapper in [WRAPPER, 'scripts/run_ESM2ne_viral.sh']:
            text = (work/wrapper).read_text()
            assert 'python -m scripts.LitESM2ne_partial_trainer' in text
            assert 'python -m scripts.LitESM2ne_Lora_trainer' not in text
            assert 'CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"' in text
            subprocess.run([shutil.which('bash'), '-n', str(work/wrapper)], check=True, capture_output=True)
        parser_scope = {'argparse': argparse,
                        'torch': SimpleNamespace(cuda=SimpleNamespace(is_available=lambda: False))}
        parse_node = next(n for n in ast.parse(historical_text).body
                          if isinstance(n, ast.FunctionDef) and n.name == 'parse_args')
        exec(compile(ast.Module(body=[parse_node], type_ignores=[]), 'historical_argument_parser', 'exec'), parser_scope)
        saved_argv = sys.argv
        try:
            sys.argv = ['restored_trainer', '-i', 'fixture.csv', '-o', 'results',
                        '--checkpoint', 'esm2_t33_650M_UR50D',
                        '--split_strategy', 'site_split', '--seed', '42']
            parsed = parser_scope['parse_args']()
            assert parsed.checkpoint == 'esm2_t33_650M_UR50D' and parsed.split_strategy == 'site_split'
        finally:
            sys.argv = saved_argv
        restored_rows = []
        for replicate in [1, 2, 3]:
            seed = lasso['make_seed']('IF1_ECOLI', replicate)
            obj = SimpleNamespace(split_strategy='site_split')
            setup_function(historical_text, seed)(obj)
            _, expected = lasso['split_data'](data, seed)
            assert set(obj.val_dataset.ID) == set(expected.ID)
            assert set(obj.train_dataset.site).isdisjoint(obj.val_dataset.site)
            restored_rows.append({'replicate': replicate, 'validation_rows': len(expected),
                                  'matches_reference_fold': True})
    return {'status': 'passed', 'meaning': 'Expected defects reproduced; local corrections and scoped historical reconstruction checks pass',
            'commit': inventory['commit'], 'integrity_files_checked': file_count,
            'original_dispatch': original_dispatch, 'patched_dispatch': fixed_dispatch,
            'original_setup_error': message, 'patched_setup_checks': rows,
            'negative_control': 'Original Lasso and LoRA helper splits agree; unfair-fold allegation is not supported',
            'missing_entrypoint': {'path': missing, 'absent_in_pinned_tree': True,
                                   'resolved_by_local_fixes': False,
                                   'restored_by_reconstruction_patch': True},
            'historical_reconstruction': {'source_commit': restoration['restoration_commit'],
                                          'historical_source_hashes_verified': True,
                                          'documented_device_selection_adaptation': 'logical device 0; visible-GPU override preserved',
                                          'wrapper_syntax_checked': True,
                                          'checkpoint_CLI_compatible': True,
                                          'training_method_replaced_by_LoRA': False,
                                          'data_setup_checks': restored_rows},
            'full_training_run': False, 'paper_results_reproduced': False,
            'external_human_confirmation': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=REPO/'artifacts/paper-code-auditor/verification.json')
    args = parser.parse_args()
    if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'):
        parser.error('Assertions must be enabled; do not use -O or PYTHONOPTIMIZE')
    result = check()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
