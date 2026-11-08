#!/usr/bin/env python3
"""Classify FASTQ reads using MC-HYD and kaijux (MIT License)."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
from itertools import zip_longest
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
INDEX = ROOT / 'kaiju/kaiju_hydrogenases_20261025.fmi'
MAPPING = ROOT / 'kaiju/protein_to_subgroup.tsv'


def fastq(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence, plus, quality = (handle.readline().rstrip('\r\n') for _ in range(3))
            fields = header[1:].split()
            if (not header.startswith('@') or not fields or not sequence
                    or not plus.startswith('+') or len(sequence) != len(quality)
                    or any(c.upper() not in 'ACGTRYSWKMBDHVN' for c in sequence)):
                raise ValueError(f'Invalid or truncated four-line FASTQ record in {path}')
            yield fields[0], sequence, quality


def pair_id(identifier):
    return re.sub(r'/[12]$', '', identifier)


def load_mapping(path):
    mapping = {}
    with path.open() as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if len(row) != 3 or row[0] in mapping or row[1] not in {'[NiFe]', '[FeFe]', '[Fe]'}:
                raise ValueError('Invalid or duplicate protein-to-subgroup mapping')
            mapping[row[0]] = (row[1], row[2])
    if not mapping:
        raise ValueError('Empty protein-to-subgroup mapping')
    return mapping


def parse_raw(path, mapping, expected, calls):
    """kaijux columns: C/U, read identifier, score, comma-separated proteins."""
    observed = set()
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2 or fields[0] not in {'C', 'U'}:
                raise ValueError('Malformed kaijux output')
            identifier = fields[1]
            if identifier not in expected or identifier in observed:
                raise ValueError(f'Unexpected or duplicate kaijux read: {identifier}')
            observed.add(identifier)
            if fields[0] == 'U':
                continue
            if len(fields) < 4:
                raise ValueError('Classified kaijux record lacks score/protein IDs')
            score = int(fields[2])
            refs = [r for r in fields[3].split(',') if r]
            if not refs:
                raise ValueError('Classified kaijux record has no reference proteins')
            unit = int(identifier.split('m')[0][1:])
            for ref in refs:
                if ref not in mapping:
                    raise ValueError(f'Kaiju reference absent from bundled mapping: {ref}')
                group = mapping[ref]
                entry = calls[unit].setdefault(group, [score, set()])
                entry[0] = max(entry[0], score)
                entry[1].add(ref)
    if observed != expected:
        raise ValueError('Kaijux did not report every input read; refusing incomplete results')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-1', '--reads1', required=True, type=Path)
    parser.add_argument('-2', '--reads2', type=Path)
    parser.add_argument('-o', '--output', required=True, type=Path)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    if not shutil.which('kaijux'):
        parser.error('kaijux is not on PATH; install Kaiju (tested with 1.10.3)')
    if not INDEX.is_file():
        parser.error(f'Missing FM-index: {INDEX}')
    mapping = load_mapping(MAPPING)
    args.output.mkdir(parents=True, exist_ok=False)
    units, names, calls = [], set(), defaultdict(dict)
    with tempfile.TemporaryDirectory(prefix='mc_hyd_reads_') as tmp:
        # Internal IDs prevent kaijux stripping slash-containing original identifiers.
        prepared = Path(tmp) / 'reads.fastq'
        with prepared.open('w') as dest, (args.output / 'read_id_map.tsv').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['kaijux_id', 'read_id', 'source_read_id', 'mate'])
            records = zip_longest(fastq(args.reads1), fastq(args.reads2)) if args.reads2 else ((r,) for r in fastq(args.reads1))
            for number, mates in enumerate(records):
                if any(r is None for r in mates):
                    raise ValueError('Paired FASTQs have different numbers of reads')
                identifier = pair_id(mates[0][0]) if args.reads2 else mates[0][0]
                if args.reads2 and pair_id(mates[1][0]) != identifier:
                    raise ValueError(f'Mismatched mate identifiers: {mates[0][0]}, {mates[1][0]}')
                if identifier in names:
                    raise ValueError(f'Duplicate input identifier: {identifier}; supply mates with -1 and -2')
                names.add(identifier)
                units.append((identifier, ','.join(r[0] for r in mates)))
                for mate, (original, sequence, quality) in enumerate(mates, 1):
                    internal = f'r{number}m{mate}'
                    dest.write(f'@{internal}\n{sequence}\n+\n{quality}\n')
                    writer.writerow([internal, identifier, original, mate])
        if not units:
            raise ValueError('FASTQ input contains no reads')
        raw = args.output / 'kaiju_raw.tsv'
        with (args.output / 'kaiju.log').open('w') as log:
            subprocess.run(['kaijux', '-f', str(INDEX), '-i', str(prepared),
                            '-a', 'greedy', '-e', '3', '-z', str(args.threads), '-o', str(raw)],
                           check=True, stdout=log, stderr=subprocess.STDOUT)
        expected = {f'r{i}m{mate}' for i in range(len(units)) for mate in range(1, 3 if args.reads2 else 2)}
        parse_raw(raw, mapping, expected, calls)
    counts = Counter()
    families = defaultdict(set)
    with (args.output / 'hydrogenase_reads.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['read_id', 'family', 'subgroup', 'score', 'protein_ids', 'source_read_ids'])
        for unit, groups in sorted(calls.items()):
            for (fam, group), (score, refs) in sorted(groups.items()):
                writer.writerow([units[unit][0], fam, group, score, ','.join(sorted(refs)), units[unit][1]])
                counts[fam, group] += 1
                families[fam].add(unit)
    with (args.output / 'hydrogenase_read_summary.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['family', 'subgroup', 'count'])
        writer.writerows((*key, count) for key, count in sorted(counts.items()))
    unit_name = 'Read pairs' if args.reads2 else 'Reads'
    print(f'MC-HYD annotation completed\n\nInput: {args.reads1}\n{unit_name} analyzed: {len(units)}'
          f'\nClassified: {len(calls)}\nUnclassified: {len(units) - len(calls)}'
          f'\nHydrogenase-associated {unit_name.lower()}: {len(calls)}')
    print('\nFamily\tCount (unique units)')
    for fam, members in sorted(families.items()):
        print(fam, len(members), sep='\t')
    print('\nFamily\tSubgroup\tCount')
    for key, count in sorted(counts.items()):
        print(*key, count, sep='\t')
    print('Ambiguous units contribute once to each supported subgroup; subgroup counts may exceed classified units.')
    print(f'\nResults saved to: {args.output.resolve()}/')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, EOFError, subprocess.CalledProcessError) as exc:
        print(f'Error: {exc}. See kaiju.log if available. Use a new output directory for each run.', file=sys.stderr)
        sys.exit(1)
