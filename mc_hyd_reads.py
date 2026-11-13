#!/usr/bin/env python3
"""MC-HYD read classification using KaijuX."""

import argparse
import csv
import os
import re
import shutil
import sqlite3
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INDEX = ROOT / "kaiju/kaiju_hydrogenases_20261025.fmi"
MAPPING = ROOT / "kaiju/protein_to_subgroup.tsv"

MATE = re.compile(r"(?<![A-Za-z0-9])R([12])(?![A-Za-z0-9])", re.I)


def available_threads():
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    if allocated:
        return max(1, allocated - 2)
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    return max(1, cpus - 2)


def find_pair(folder):
    pairs = {}

    for path in folder.iterdir():
        if not path.is_file() or not path.name.lower().endswith(".fastq.gz"):
            continue

        matches = list(MATE.finditer(path.name))
        if not matches:
            continue
        if len(matches) != 1:
            raise ValueError(f"Ambiguous mate name: {path.name}")

        match = matches[0]
        key = path.name[:match.start()] + "{mate}" + path.name[match.end():]
        mate = int(match.group(1))
        pairs.setdefault(key, {})
        if mate in pairs[key]:
            raise ValueError(f"Duplicate R{mate} for {key}")
        pairs[key][mate] = path

    if len(pairs) != 1 or set(next(iter(pairs.values()), {})) != {1, 2}:
        raise ValueError(f"Expected exactly one R1/R2 FASTQ.gz pair in {folder}")

    pair = next(iter(pairs.values()))
    return pair[1], pair[2]


def load_mapping():
    with MAPPING.open() as handle:
        return {
            protein: (family, subgroup)
            for protein, family, subgroup in csv.reader(handle, delimiter="\t")
        }


def classify(fastq, threads, mapping, db, log):
    cmd = [
        "kaijux", "-f", str(INDEX), "-i", str(fastq),
        "-a", "greedy", "-e", "3", "-z", str(threads)
    ]

    total = 0

    with log.open("w") as err:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=err,
            text=True, bufsize=1024 * 1024
        )

        try:
            with db:
                for line in proc.stdout:
                    total += 1
                    fields = line.rstrip("\n").split("\t")

                    if len(fields) < 2 or fields[0] not in ("C", "U"):
                        raise ValueError(f"Invalid KaijuX output at record {total}")

                    if fields[0] == "U":
                        continue

                    if len(fields) < 4:
                        raise ValueError(f"Incomplete KaijuX hit at record {total}")

                    read_id = re.sub(r"/[12]$", "", fields[1])
                    score = int(fields[2])
                    proteins = [p for p in fields[3].split(",") if p]

                    if not proteins:
                        raise ValueError(f"Missing reference IDs at record {total}")

                    groups = {}
                    for protein in proteins:
                        if protein not in mapping:
                            raise ValueError(f"Unknown reference: {protein}")
                        groups.setdefault(mapping[protein], set()).add(protein)

                    for (family, subgroup), refs in groups.items():
                        db.execute(
                            """
                            INSERT INTO calls VALUES (?, ?, ?, ?, ?)
                            ON CONFLICT(read_id, family, subgroup)
                            DO UPDATE SET
                                score = MAX(score, excluded.score),
                                refs = refs || ',' || excluded.refs
                            """,
                            (read_id, family, subgroup, score,
                             ",".join(sorted(refs)))
                        )

            if proc.wait():
                raise RuntimeError(f"KaijuX failed on {fastq}; see {log}")

        except BaseException:
            if proc.poll() is None:
                proc.kill()
            proc.wait()
            raise

        finally:
            proc.stdout.close()

    return total


def run_sample(name, r1, r2, output, threads, mapping):
    destination = output / name
    destination.mkdir()

    scratch = os.environ.get("SLURM_TMPDIR") or None

    with tempfile.TemporaryDirectory(prefix="mc_hyd_", dir=scratch) as tmp:
        with sqlite3.connect(Path(tmp) / "hits.sqlite") as db:
            db.execute("PRAGMA journal_mode=OFF")
            db.execute("PRAGMA synchronous=OFF")
            db.execute("""
                CREATE TABLE calls (
                    read_id TEXT,
                    family TEXT,
                    subgroup TEXT,
                    score INTEGER,
                    refs TEXT,
                    PRIMARY KEY (read_id, family, subgroup)
                )
            """)

            print(f"{name}: R1", flush=True)
            n1 = classify(
                r1, threads, mapping, db, destination / "R1.log"
            )

            print(f"{name}: R2", flush=True)
            n2 = classify(
                r2, threads, mapping, db, destination / "R2.log"
            )

            if n1 != n2:
                raise ValueError(f"{name}: R1/R2 read counts differ")

            counts = Counter()

            with (destination / "hydrogenase_reads.tsv").open(
                "w", newline=""
            ) as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow([
                    "read_id", "family", "subgroup", "score", "protein_ids"
                ])

                for read_id, family, subgroup, score, refs in db.execute(
                    "SELECT * FROM calls ORDER BY read_id, family, subgroup"
                ):
                    writer.writerow([
                        read_id, family, subgroup, score,
                        ",".join(sorted(set(refs.split(","))))
                    ])
                    counts[family, subgroup] += 1

            positive = db.execute(
                "SELECT COUNT(DISTINCT read_id) FROM calls"
            ).fetchone()[0]

    with (destination / "hydrogenase_read_summary.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["family", "subgroup", "count"])
        for (family, subgroup), count in sorted(counts.items()):
            writer.writerow([family, subgroup, count])

    print(f"{name}: {n1:,} pairs; {positive:,} positive", flush=True)
    return counts, n1, positive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=available_threads())
    args = parser.parse_args()

    if args.threads < 1:
        parser.error("--threads must be positive")
    if not shutil.which("kaijux"):
        parser.error("kaijux not found")
    if not INDEX.is_file() or not MAPPING.is_file():
        parser.error("MC-HYD database files missing")
    if not args.input.is_dir():
        parser.error("Input folder does not exist")
    if args.output.exists():
        parser.error("Output folder already exists")

    mapping = load_mapping()

    fastqs = list(args.input.glob("*.fastq.gz"))

    if fastqs:
        samples = [(args.input.name, *find_pair(args.input))]
    else:
        samples = [
            (folder.name, *find_pair(folder))
            for folder in sorted(args.input.iterdir())
            if folder.is_dir()
            and any(folder.glob("*.fastq.gz"))
        ]

    if not samples:
        parser.error("No paired FASTQ.gz files found")

    args.output.mkdir(parents=True)

    results = {}
    totals = {}

    for name, r1, r2 in samples:
        counts, total, positive = run_sample(
            name, r1, r2, args.output, args.threads, mapping
        )
        results[name] = counts
        totals[name] = (total, positive)

    with (args.output / "hydrogenase_type_by_sample.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["family", "subgroup", *results])

        for family, subgroup in sorted(set(mapping.values())):
            writer.writerow([
                family, subgroup,
                *(results[s].get((family, subgroup), 0) for s in results)
            ])

    with (args.output / "sample_summary.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "sample", "read_pairs", "hydrogenase_positive_pairs"
        ])
        for name, (total, positive) in totals.items():
            writer.writerow([name, total, positive])

    print(f"\nCompleted {len(results)} samples. Results: {args.output}")


if __name__ == "__main__":
    main()
    
