# MC-HYD: a curated hydrogenase database and annotation toolkit

[![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/powered-by-coffee.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/60-percent-of-the-time-works-every-time.svg)](https://forthebadge.com)

[![giovannellilab](https://img.shields.io/badge/BY-Giovannelli_Lab-blue)](https://www.donatogiovannelli.com)
[![funded-by-erc](https://img.shields.io/badge/Funded%20by-ERC-ff6400.svg)](https://erc.europa.eu/homepage)
[![project-coevolve](https://img.shields.io/badge/Project-ERC%20CoEvolve-000fa9.svg)](https://www.coevolve.eu/)


MC-HYD is the Giovannelli Lab's shared, versioned resource for consistent
hydrogenase annotation across research projects. Database release **20261025**
contains **7,134 curated hydrogenase catalytic-subunit proteins**, **39
subgroup-specific HMMs**, and a Kaiju FM-index built from **7,148 normalized
protein segments**. The difference between protein and segment counts is
intentional: normalization can produce multiple indexed segments per protein.

Two independent commands identify and classify hydrogenases in the three
families **[NiFe]**, **[FeFe]**, and **[Fe]**:

- `mc_hyd_annotate.py`: HMMER annotation of protein FASTA files, including
  predicted proteomes from genomes and metagenome-assembled genomes (MAGs).
- `mc_hyd_reads.py`: Kaiju annotation directly from metagenomic FASTQ reads.
  Assignments describe hydrogenase functional subgroups, not organism taxonomy.

| Family | Subgroups | HMMs |
|---|---|---:|
| [NiFe] | 1a–1l, 2a–2e, 3a–3d, 4a–4i | 30 |
| [FeFe] | A1–A4, B, C1–C3 | 8 |
| [Fe] | Fe | 1 |

## Installation

Clone the repository; both Python scripts use only Python's standard library
(Python 3.8 or later). Install the external tool for the workflow you need.
HMMER and Kaiju need not both be installed.

```bash
git clone https://github.com/giovannellilab/MC_HYD.git
cd MC_HYD
python --version
```

### Protein annotation: HMMER 3

Ubuntu/Debian:

```bash
sudo apt update
sudo apt install python3 hmmer
```

Or Conda/Bioconda:

```bash
conda create -n mc-hyd-hmm -c conda-forge -c bioconda python=3.11 hmmer
conda activate mc-hyd-hmm
```

Verify the executables:

```bash
command -v hmmscan
command -v hmmpress
hmmscan -h
```

See the [official HMMER website](https://www.hmmer.org/),
[installation instructions](https://www.hmmer.org/install.html), and
[documentation](https://www.hmmer.org/documentation.html).

### Read annotation: Kaiju

```bash
conda create -n mc-hyd-reads -c conda-forge -c bioconda python=3.11 kaiju=1.10.3
conda activate mc-hyd-reads
command -v kaijux
kaijux -h
```

MC-HYD requires **`kaijux`**, not the conventional `kaiju` taxonomy command.
Kaiju 1.10.3 was used for the local smoke tests. The
[Bioconda build recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/kaiju/build.sh)
installs the `kaiju*` binaries, including `kaijux`; check its availability in
your activated environment as shown above. For source installation and full
usage documentation, see the
[official Kaiju repository](https://github.com/bioinformatics-centre/kaiju#installation)
and [usage instructions](https://github.com/bioinformatics-centre/kaiju#usage).
No taxonomy dumps or additional database downloads are needed for MC-HYD.

## Protein FASTA annotation

```bash
python mc_hyd_annotate.py proteins.faa -o results/MAG_001 --threads 8
```

Input must be an amino-acid FASTA with nonempty sequences and unique identifiers
(the first whitespace-delimited token after `>`). Original identifiers are
preserved; descriptions are not part of the identifier. Selenocysteine (`U`) is
converted to cysteine (`C`) in a temporary query file, consistent with model
development. The input file is unchanged.

The script prepares HMMER's binary indexes with `hmmpress` when missing or older
than the bundled library, then runs `hmmscan` with a **400-bit full-sequence
score threshold** (`-T 400`). It retains the highest-scoring model at or above
400 bits for each protein and reports that model's subgroup and family.
The table's E-value is the full-sequence E-value. Equal-score ties retain the
first hit in the HMMER table. There are no additional sequence filters or score
normalizations. The library directory must be writable for initial indexing;
prepare it once before simultaneous runs on shared installations.

Outputs in `results/MAG_001/`:

| File | Contents |
|---|---|
| `hydrogenases.tsv` | `protein_id`, `family`, `subgroup`, `bitscore`, `evalue` |
| `hydrogenase_summary.tsv` | `family`, `subgroup`, `count` |
| `hmmscan.tbl` | Original HMMER tabular results |
| `hmmscan.log` | HMMER report and diagnostics |

Each accepted protein contributes exactly once to the summary. A protein absent
from the predictions has no accepted hit under this operational threshold;
absence is not proof of biological absence in an incomplete genome or proteome.

## Metagenomic read annotation

Single-end (plain FASTQ or gzip-compressed FASTQ):

```bash
python mc_hyd_reads.py -1 sample.fastq.gz -o results/sample_001 --threads 8
```

Paired-end:

```bash
python mc_hyd_reads.py -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    -o results/sample_001_paired --threads 8
```

Input uses conventional four-line FASTQ records. Paired files must have the
same number of records in matching order, with matching first-token identifiers
after removal of terminal `/1` or `/2`. Illumina mate information after whitespace
is accepted. Identifiers must be unique within the counting unit. Supply mates
using `-1` and `-2`; do not supply an interleaved paired file as single-end input.

The original notebook command, reused here, is:

```text
kaijux -f kaiju/kaiju_hydrogenases_20261025.fmi -i READS \
    -a greedy -e 3 -z THREADS -o OUTPUT
```

Other search settings retain the Kaiju 1.10.3 defaults: minimum fragment length
11, minimum greedy score 65, E-value cutoff 0.01, and SEG low-complexity filtering
enabled. Use the pinned version above to reproduce these defaults. The HMMER
400-bit threshold does **not** apply to Kaiju scores.

The established benchmark searches mates individually and takes the union of
supported subgroups per pair. MC-HYD preserves this procedure by preparing one
temporary single-end stream containing both mates, rather than switching to
Kaiju's native joint paired search. For `-1` only, the counting unit is one read;
for `-1` plus `-2`, it is one pair. Both mates supporting the same subgroup count
once. Multiple protein hits in a subgroup also count once. If different
subgroups are supported, every supported subgroup is retained without an
arbitrary tie-break. No abundance normalization is performed.

`kaijux` output is parsed as classification flag, read ID, score, and a
comma-separated list of reference protein IDs. The bundled
`protein_to_subgroup.tsv` maps those IDs to functional family and subgroup.
Unknown reference IDs or incomplete raw output stop the run rather than silently
losing assignments.

| File | Contents |
|---|---|
| `hydrogenase_reads.tsv` | `read_id`, `family`, `subgroup`, `score`, `protein_ids`, `source_read_ids` |
| `hydrogenase_read_summary.tsv` | `family`, `subgroup`, `count` |
| `kaiju_raw.tsv` | Original, unmodified `kaijux` output |
| `read_id_map.tsv` | Internal `kaijux_id`, output `read_id`, original `source_read_id`, and `mate` |
| `kaiju.log` | External-tool diagnostics |

There is one prediction row per supported read/subgroup or pair/subgroup.
`score` is the highest observed Kaiju score supporting that subgroup across the
unit; `protein_ids` lists its supporting references. Scores are Kaiju alignment
scores, not HMMER bit scores. `source_read_ids` preserves original first-token
read identifiers; paired output `read_id` removes terminal mate suffixes.
Internal IDs prevent Kaiju's own slash-based name truncation from merging
unrelated reads. `read_id_map.tsv` connects raw output to the original IDs.

The terminal reports analyzed, classified, unclassified, and hydrogenase-associated
units. Classified units are those with any mapped hydrogenase call. Family counts
are unique units per family. Subgroup counts equal the rows in the prediction
table for each subgroup. Ambiguous units contribute to multiple subgroups (and
potentially families), so those counts can sum to more than the classified total.
Read calls identify hydrogenase-associated sequence; they do not establish a
complete catalytic gene or conventional taxonomic assignment.

## Running and troubleshooting

Both scripts resolve bundled references relative to their own location and work
from any current directory when invoked by absolute path. `--threads` must be
positive. Use `--help` for the command-line interface. Each run requires a **new
output directory** to prevent stale results being mistaken for a completed run.
Parent directories are created automatically. Failed runs return a nonzero exit
status and may leave raw files or logs for inspection; only trust runs reporting
`MC-HYD annotation completed`. Runs with no accepted hits produce header-only
prediction and summary tables. Temporary query files are removed automatically.
The read wrapper holds identifiers and accepted assignments in memory; allow
memory proportional to the number of input units for large metagenomes.

## Validation and limitations

Established HMM validation at the 400-bit threshold found:

- 99.86% sensitivity on internal hydrogenase references.
- 99.73% sensitivity on independent external hydrogenases.
- Rejection of all 422 curated confounders.
- Rejection of all 12,482 non-target GS-24-all proteins.
- In an annotation-based screen of 576,073 reviewed UniProtKB/Swiss-Prot proteins,
  59 were accepted and no unrelated-protein false positives were identified.

These are validation-set results, not a guarantee for every sequence or sample.
The Swiss-Prot screen depends on available annotations and their interpretation;
it is not an exhaustive experimentally labelled specificity experiment.

Historical read-based results were verified from the development project's
`results/external_art/{read_pair,protein}_classification.tsv` and
`results/confounder_art/{read_pair,protein}_classification.tsv`, filtering
`method == kaiju`:

| Historical Kaiju test | Correct only | Ambiguous | Wrong only | Undetected |
|---|---:|---:|---:|---:|
| External simulated read pairs (56,010) | 55,495 | 220 | 295 | 0 |
| External source proteins (372) | 333 | 37 | 2 | 0 |

An ambiguous call includes the expected subgroup and at least one additional
subgroup; expected external `C` annotations accepted C1–C3 as in the original
analysis. Protein-level calls combine all simulated pairs from each protein.
In the curated confounder simulation, 2,997 of 54,015 pairs had false calls;
51,018 were rejected. At the source-protein level, 105 of 422 had at least one
false call and 317 were rejected. Read-based confounder rejection is therefore
not equivalent to HMM-level rejection.

These are historical benchmark results, not a new validation of the final
7,134-protein/7,148-segment index: the development notebooks used an earlier
reference build. The final reference assets were supplied as finalized and were
copied without alteration or rebuilding. Release smoke tests exercise the final
assets and wrappers on small existing-data subsets; they do not constitute a
repeat benchmark or establish metagenomic abundance accuracy.

## Repository contents

```text
MC_HYD/
├── README.md
├── LICENSE
├── LICENSE-DATA.md
├── .gitignore
├── mc_hyd_annotate.py
├── mc_hyd_reads.py
├── database/
│   └── MC_HYD_DB_20261025.fasta
├── hmm/
│   └── MC_HYD_DB_HMM_20261025.hmm
└── kaiju/
    ├── kaiju_hydrogenases_20261025.fmi
    └── protein_to_subgroup.tsv
```

All four reference assets are tracked directly in Git. Generated HMMER indexes,
analysis results, environments, and development notebooks are not release assets.

## Provenance and citation

MC-HYD builds on the original HydDB classification framework, with additional
sequences from NCBI RefSeq/NR and other publicly available sources. Original
reference accessions and source annotations are retained in the reference FASTA.
The MC-HYD database and validation manuscript is **in preparation**.

Please cite MC-HYD by repository URL and release, together with the relevant
foundational resources:

- Søndergaard, D., Pedersen, C. N. S., & Greening, C. (2016). HydDB: A web tool
  for hydrogenase classification and analysis. *Scientific Reports*, 6, 34212.
  https://doi.org/10.1038/srep34212
- Eddy, S. R. (2011). Accelerated Profile HMM Searches. *PLoS Computational
  Biology*, 7(10), e1002195. https://doi.org/10.1371/journal.pcbi.1002195
- Menzel, P., Ng, K. L., & Krogh, A. (2016). Fast and sensitive taxonomic
  classification for metagenomics with Kaiju. *Nature Communications*, 7, 11257.
  https://doi.org/10.1038/ncomms11257

In publications, report the MC-HYD release and database version, annotation
method and tool version, and relevant threshold or parameters. For example:
“Proteins were annotated using MC-HYD database 20261025 and HMMER hmmscan,
retaining the highest-scoring subgroup per protein at ≥400 full-sequence bits.”
For reads, report `kaijux -a greedy -e 3`, the Kaiju version, single-read or
paired counting, and retention of ambiguous subgroup assignments.

## Licenses

Original Python scripts are distributed under the standard **MIT License** in
[LICENSE](LICENSE). Original MC-HYD curation, annotations, organization, and
derived models are licensed under **CC BY 4.0**, to the extent applicable; the
standard text and scope are in [LICENSE-DATA.md](LICENSE-DATA.md). Attribute
MC-HYD and upstream sources. Individual amino-acid sequences remain subject to
applicable upstream rights; MC-HYD does not claim exclusive ownership of them.
HMMER and Kaiju are separately licensed third-party software, installed by the
user and not bundled with this repository.
