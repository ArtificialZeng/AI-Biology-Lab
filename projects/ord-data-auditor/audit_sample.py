#!/usr/bin/env python3
"""Audit only CAS fields in a frozen, hash-checked ORD convenience sample."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import urllib.request

from cas_guard import inspect_cas

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]


def walk(message, path=''):
    yield path, message
    for field, value in message.ListFields():
        if field.type != field.TYPE_MESSAGE:
            continue
        child = path+'.'+field.name if path else field.name
        if field.label == field.LABEL_REPEATED:
            if field.message_type.GetOptions().map_entry:
                if field.message_type.fields_by_name['value'].type == field.TYPE_MESSAGE:
                    for key, item in value.items():
                        yield from walk(item, f'{child}[{key}]')
            else:
                for index, item in enumerate(value):
                    yield from walk(item, f'{child}[{index}]')
        else:
            yield from walk(value, child)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=REPO/'artifacts/ord-data-auditor/datasets')
    parser.add_argument('--download', action='store_true', help='Fetch the frozen 2.4MB sample if absent')
    parser.add_argument('--output', type=Path, default=REPO/'artifacts/ord-data-auditor/sample_audit.json')
    args = parser.parse_args()
    import fastparquet
    from ord_schema.proto import reaction_pb2
    manifest = json.loads((ROOT/'evidence/sample_manifest.json').read_text())
    args.data_dir.mkdir(parents=True, exist_ok=True)

    def obtain(entry):
        path = args.data_dir/entry['filename']
        if not path.exists():
            if not args.download:
                raise FileNotFoundError(f'{entry["filename"]} absent; use --download or --data-dir')
            with urllib.request.urlopen(entry['url'], timeout=60) as response:
                content = response.read(entry['bytes']+1)
            if len(content) != entry['bytes'] or hashlib.sha256(content).hexdigest() != entry['sha256']:
                raise RuntimeError(f'Download hash/size mismatch: {entry["filename"]}')
            path.write_bytes(content)
        if path.stat().st_size != entry['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            raise RuntimeError(f'Input hash/size mismatch: {entry["filename"]}')
        return path

    with ThreadPoolExecutor(max_workers=4) as pool:
        paths = list(pool.map(obtain, manifest['files']))
    retired_file = ROOT/'fixtures/data/retired_datasets.csv'
    if hashlib.sha256(retired_file.read_bytes()).hexdigest() != manifest['retired_mapping_sha256']:
        raise RuntimeError('Frozen retirement mapping hash mismatch')
    retired = {r['retired_dataset_id'] for r in csv.DictReader(retired_file.open())}
    counts = Counter(); findings = []; datasets = []
    for path in sorted(paths):
        table = fastparquet.ParquetFile(path)
        dataset_id = table.key_value_metadata.get('ord.dataset_id', '')
        if dataset_id in retired:
            raise RuntimeError(f'Frozen sample unexpectedly includes a retired ID: {dataset_id}')
        frame = table.to_pandas()
        if list(frame.columns) != ['reaction_id', 'reaction']:
            raise RuntimeError(f'Unexpected table columns: {path.name}')
        counts['datasets'] += 1
        datasets.append({'filename': path.name, 'dataset_id': dataset_id, 'rows': len(frame)})
        for row in frame.itertuples(index=False):
            reaction = reaction_pb2.Reaction.FromString(row.reaction)
            if row.reaction_id != reaction.reaction_id:
                raise RuntimeError(f'Row/Protobuf reaction ID mismatch: {path.name}')
            counts['reactions'] += 1
            for position, message in walk(reaction):
                if message.DESCRIPTOR.name != 'CompoundIdentifier' or message.type != message.CAS_NUMBER:
                    continue
                counts['CAS_identifiers'] += 1
                check = inspect_cas(message.value)
                if check['status'] != 'checksum_valid':
                    counts[check['status']] += 1
                    findings.append(dict(check, dataset_file=path.name,
                                         reaction_id=row.reaction_id, path=position))
    observed = {key: counts[key] for key in manifest['expected_counts']}
    if observed != manifest['expected_counts']:
        raise RuntimeError(f'Frozen-sample reference mismatch: {observed}')
    report = {'status': 'passed', 'scope': 'CAS-only audit; not complete ORD semantic or reaction validation',
              'mirror_commit': manifest['mirror_commit'], 'sample_selection': manifest['selection'],
              'counts': observed, 'datasets': datasets, 'findings': findings,
              'reader': {'name': 'fastparquet', 'version': fastparquet.__version__},
              'known_existing_format_warnings': 9, 'checksum_gap_occurrences': 1,
              'automatic_data_changes': 0, 'registry_identity_verified': False,
              'population_error_rate_estimated': False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in {'datasets','findings'}}, indent=2))


if __name__ == '__main__':
    main()
