#!/usr/bin/env python3
"""Annotate protein FASTA with the MC-HYD HMM library (MIT License)."""
import argparse
from collections import Counter
import csv
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
HMM = ROOT / 'hmm/MC_HYD_DB_HMM_20261025.hmm'


def family(group):
    if group == 'Fe':
        return '[Fe]'
    if group in {'A1', 'A2', 'A3', 'A4', 'B', 'C1', 'C2', 'C3'}:
        return '[FeFe]'
    if group in {f'{n}{c}' for n, letters in [(1, 'abcdefghijkl'), (2, 'abcde'),
                                             (3, 'abcd'), (4, 'abcdefghi')]
                 for c in letters}:
        return '[NiFe]'
    raise ValueError(f'Unknown hydrogenase subgroup: {group}')


def parse_hits(path, identifiers):
    """Use full-sequence E-value and score, not domain scores."""
    hits = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 18:
                raise ValueError('Malformed hmmscan table')
            group, protein = fields[0], fields[2]
            score, evalue = float(fields[5]), float(fields[4])
            if protein not in identifiers:
                raise ValueError(f'Unexpected HMMER query: {protein}')
            fam = family(group)
            if score >= 400 and (protein not in hits or score > hits[protein][2]):
                hits[protein] = (fam, group, score, evalue)
    return hits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('proteins', type=Path)
    parser.add_argument('-o', '--output', required=True, type=Path)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    for exe in ('hmmscan', 'hmmpress'):
        if not shutil.which(exe):
            parser.error(f'{exe} is not on PATH; install HMMER 3')
    if not HMM.is_file():
        parser.error(f'Missing HMM library: {HMM}')
    # A new directory prevents stale result tables from masquerading as a successful run.
    args.output.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='mc_hyd_') as tmp:
        normalized = Path(tmp) / 'proteins.faa'
        identifiers = []
        seen = set()
        has_sequence = False
        with args.proteins.open() as source, normalized.open('w') as dest:
            for line in source:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if identifiers and not has_sequence:
                        raise ValueError(f'Empty sequence: {identifiers[-1]}')
                    fields = line[1:].split()
                    if not fields or fields[0] in seen:
                        raise ValueError('Empty or duplicate FASTA identifier')
                    identifier = fields[0]
                    identifiers.append(identifier)
                    seen.add(identifier)
                    has_sequence = False
                    dest.write(f'>{identifier}\n')
                else:
                    if not identifiers or not all(c in 'ACDEFGHIKLMNPQRSTVWYUBZXJO*' for c in line.upper()):
                        raise ValueError('Invalid amino-acid FASTA sequence')
                    dest.write(line.upper().replace('U', 'C') + '\n')
                    has_sequence = True
        if not identifiers or not has_sequence:
            raise ValueError('FASTA contains no proteins or ends with an empty sequence')
        indexes = [Path(str(HMM) + suffix) for suffix in ('.h3f', '.h3i', '.h3m', '.h3p')]
        if any(not p.exists() or p.stat().st_mtime < HMM.stat().st_mtime for p in indexes):
            subprocess.run(['hmmpress', '-f', str(HMM)], check=True, stdout=subprocess.DEVNULL)
        raw = args.output / 'hmmscan.tbl'
        with (args.output / 'hmmscan.log').open('w') as log:
            subprocess.run(['hmmscan', '--cpu', str(args.threads), '--noali', '-T', '400',
                            '--tblout', str(raw), str(HMM), str(normalized)],
                           check=True, stdout=log, stderr=subprocess.STDOUT)
        hits = parse_hits(raw, seen)
        counts = Counter((v[0], v[1]) for v in hits.values())
        with (args.output / 'hydrogenases.tsv').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['protein_id', 'family', 'subgroup', 'bitscore', 'evalue'])
            for protein in identifiers:
                if protein in hits:
                    writer.writerow([protein, *hits[protein]])
        with (args.output / 'hydrogenase_summary.tsv').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['family', 'subgroup', 'count'])
            writer.writerows((*key, count) for key, count in sorted(counts.items()))
    print(f'MC-HYD annotation completed\n\nInput: {args.proteins}\nProteins analyzed: {len(identifiers)}'
          f'\nHydrogenases detected: {len(hits)}')
    print('\nFamily\tCount')
    for fam, count in sorted(Counter(v[0] for v in hits.values()).items()):
        print(f'{fam}\t{count}')
    print('\nFamily\tSubgroup\tCount')
    for key, count in sorted(counts.items()):
        print(*key, count, sep='\t')
    print(f'\nResults saved to: {args.output.resolve()}/')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f'Error: {exc}. See hmmscan.log if available. Use a new output directory for each run.', file=sys.stderr)
        sys.exit(1)
