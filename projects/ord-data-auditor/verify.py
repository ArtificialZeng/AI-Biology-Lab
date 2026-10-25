#!/usr/bin/env python3
"""Offline regression checks for CAS diagnostics and the proposed ORD branch."""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace

from cas_guard import inspect_cas, inspect_records

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]


def cas_branch(text):
    function = next(n for n in ast.parse(text).body
                    if isinstance(n, ast.FunctionDef) and n.name == '_validate_compound_identifier')
    node = next(n for n in ast.walk(function) if isinstance(n, ast.If)
                and isinstance(n.test, ast.Compare)
                and any(isinstance(v, ast.Attribute) and v.attr == 'CAS_NUMBER'
                        for v in n.test.comparators))
    wrapper = ast.parse('def check_branch(message, context):\n    pass\n').body[0]
    wrapper.body = node.body
    module = ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[]))
    namespace = {'re': re}
    exec(compile(module, 'pinned_CAS_branch_only', 'exec'), namespace)

    def call(value):
        warnings = []
        message = SimpleNamespace(value=value)
        namespace['check_branch'](message, SimpleNamespace(warn=warnings.append))
        assert message.value == value
        return warnings
    return call


def check():
    manifest = json.loads((ROOT/'evidence/input_manifest.json').read_text())
    for name, expected in manifest['sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, name
    valid = ['7732-18-5', '64-17-5', '58-08-2', '50-78-2']
    malformed = ['', 'unknown', '64175', '01.03.5337', '91247-08-04',
                 '1-18-5', '12345678-18-5', '7732-18-50', '7732-18-5 ', '６４-１７-５']
    mismatch = ['7732-18-4', '977739-46-3']
    for value in valid:
        assert inspect_cas(value)['status'] == 'checksum_valid'
    for value in malformed:
        assert inspect_cas(value)['status'] == 'malformed'
    for value in mismatch:
        assert inspect_cas(value)['status'] == 'checksum_mismatch'
    assert inspect_cas(None)['status'] == 'invalid_type'
    assert inspect_cas(123)['status'] == 'invalid_type'
    witness = inspect_cas('977739-46-3')
    assert witness['weighted_sum'] == 251 and witness['arithmetic_expected_digit'] == 1
    assert witness['automatic_replacement'] is None and not witness['registry_identity_verified']
    # Independent left-to-right reference arithmetic, on artificial prefixes.
    generator = random.Random(20261026)
    property_cases = []
    for length in range(2, 8):
        for _ in range(20):
            prefix = ''.join(str(generator.randrange(10)) for _ in range(length))
            second = ''.join(str(generator.randrange(10)) for _ in range(2))
            body = prefix+second
            expected = sum((len(body)-index)*int(char) for index, char in enumerate(body)) % 10
            for digit in range(10):
                value = f'{prefix}-{second}-{digit}'
                result = inspect_cas(value)
                assert result['checksum_valid'] == (digit == expected)
                assert result['automatic_replacement'] is None
                property_cases.append((value, digit == expected))
    fixtures = [json.loads(line) for line in (ROOT/'fixtures/data/cas_findings.jsonl').read_text().splitlines()]
    expected = [row['expected_status'] for row in fixtures]
    diagnosed = inspect_records(fixtures)
    assert [row['status'] for row in diagnosed] == expected
    assert Counter(expected) == {'malformed': 9, 'checksum_mismatch': 1}
    original = cas_branch((ROOT/'fixtures/upstream/ord_schema/validations.py').read_text())
    assert sum(bool(original(row['value'])) for row in fixtures) == 9
    assert original('977739-46-3') == []
    with tempfile.TemporaryDirectory(prefix='ord-cas-check-') as temporary:
        work = Path(temporary)
        shutil.copytree(ROOT/'fixtures/upstream', work, dirs_exist_ok=True)
        patch = ROOT/'patches/ord-schema-cas-checksum.patch'
        subprocess.run(['git', 'apply', '--check', str(patch)], cwd=work, check=True, capture_output=True)
        subprocess.run(['git', 'apply', str(patch)], cwd=work, check=True, capture_output=True)
        patched = cas_branch((work/'ord_schema/validations.py').read_text())
        assert sum(bool(patched(row['value'])) for row in fixtures) == 10
        assert 'invalid check digit' in patched('977739-46-3')[0]
        for value in valid:
            assert patched(value) == []
        for value in malformed+mismatch:
            if value.isascii():
                assert len(patched(value)) == 1
        # The candidate intentionally retains the upstream format expression.
        # The standalone helper is ASCII-strict; this legacy behavior is not
        # silently changed as part of an upstream checksum-only patch.
        assert patched('６４-１７-５') == original('６４-１７-５') == []
        for value, correct in property_cases:
            assert bool(patched(value)) == (not correct)
        ast.parse((work/'ord_schema/validations_test.py').read_text())
    return {'status': 'passed', 'meaning': 'Scoped identifier checks and candidate patch pass offline regression',
            'integrity_files_checked': len(manifest['sha256']),
            'known_positive_controls': len(valid), 'malformed_controls': len(malformed),
            'checksum_negative_controls': len(mismatch), 'invalid_type_controls': 2,
            'artificial_arithmetic_cases': len(property_cases),
            'fixture_counts': dict(Counter(expected)),
            'original_CAS_branch_fixture_warnings': 9, 'patched_CAS_branch_fixture_warnings': 10,
            'candidate_patch_applies': True, 'candidate_test_file_syntax_valid': True,
            'upstream_legacy_format_behavior_preserved': True,
            'input_values_modified': False, 'automatic_identity_corrections': 0,
            'full_ORD_semantic_validation_run': False, 'full_upstream_test_suite_run': False,
            'scope': 'Actual extracted CAS-specific branch only; no RDKit, reaction validation, or registry lookup'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=REPO/'artifacts/ord-data-auditor/verification.json')
    args = parser.parse_args()
    if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'):
        parser.error('Assertions must be enabled')
    result = check()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
