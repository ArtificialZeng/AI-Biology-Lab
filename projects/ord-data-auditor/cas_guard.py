#!/usr/bin/env python3
"""CAS syntax/checksum diagnostics. Never verifies or rewrites chemical identity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

FORMAT = re.compile(r'([0-9]{2,7})-([0-9]{2})-([0-9])')


def inspect_cas(value):
    result = {'value': value, 'format_valid': False, 'checksum_valid': None,
              'registry_identity_verified': False, 'automatic_replacement': None}
    if not isinstance(value, str):
        return dict(result, status='invalid_type')
    match = FORMAT.fullmatch(value)
    if match is None:
        return dict(result, status='malformed')
    body = match.group(1) + match.group(2)
    weighted_sum = sum(weight * int(digit)
                       for weight, digit in enumerate(reversed(body), start=1))
    expected = weighted_sum % 10
    valid = expected == int(match.group(3))
    result.update(format_valid=True, checksum_valid=valid,
                  status='checksum_valid' if valid else 'checksum_mismatch',
                  weighted_sum=weighted_sum, arithmetic_expected_digit=expected)
    return result


def inspect_records(records):
    output = []
    for index, record in enumerate(records, start=1):
        if record.get('identifier_type', 'CAS_NUMBER') != 'CAS_NUMBER':
            continue
        result = inspect_cas(record.get('value'))
        result.update({k: record[k] for k in ['dataset_file', 'reaction_id', 'path'] if k in record})
        result['input_row'] = index
        output.append(result)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--value', action='append', help='A value to diagnose; repeat for multiple values')
    source.add_argument('--input', type=Path, help='JSONL identifier records with value/type/context fields')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = ([{'value': value} for value in args.value] if args.value is not None else
               [json.loads(line) for line in args.input.read_text().splitlines() if line.strip()])
    results = inspect_records(records)
    report = {'scope': 'Syntax and arithmetic only; a passing checksum does not establish a registered compound',
              'checked': len(results), 'findings': sum(r['status'] != 'checksum_valid' for r in results),
              'results': results}
    text = json.dumps(report, indent=2, ensure_ascii=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end='')
    return 1 if report['findings'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
