
# MC-HYD:A curated hydrogenase database and annotation toolkit

[![forthebadge](https://forthebadge.com/api/badges/made-with-python.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/powered-by-coffee.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/api/badges/60-percent-of-the-time-works-every-time.svg)](https://forthebadge.com)

[![giovannellilab](https://img.shields.io/badge/BY-Giovannelli_Lab-blue)](https://www.donatogiovannelli.com)
[![funded-by-erc](https://img.shields.io/badge/Funded%20by-ERC-ff6400.svg)](https://erc.europa.eu/homepage)
[![project-coevolve](https://img.shields.io/badge/Project-ERC%20CoEvolve-000fa9.svg)](https://www.coevolve.eu/)


MC-HYD provides a curated reference database and validated tools for identifying and classifying hydrogenase catalytic subunits in genomes, metagenome-assembled genomes (MAGs), and metagenomic sequencing reads.

The toolkit supports two complementary approaches:

- **Profile Hidden Markov Models (HMMER)** for protein sequences, predicted ORFs, genomes, and MAGs.
- **Kaiju** for direct hydrogenase classification from metagenomic sequencing reads.

MC-HYD builds upon the original HydDB classification framework, expanding hydrogenase diversity and providing experimentally benchmarked annotation workflows.

## Database overview

The current MC-HYD release contains **7,134 curated hydrogenase catalytic-subunit sequences**, representing 39 subgroups across the three major hydrogenase families.

| Hydrogenase family | Subgroups | Models |
|---|---|---:|
| [NiFe] | 1a–1l, 2a–2e, 3a–3d, 4a–4i | 30 |
| [FeFe] | A1–A4, B, C1–C3 | 8 |
| [Fe] | Fe | 1 |
| **Total** | | **39** |

The database incorporates sequences from the original HydDB reference, NCBI RefSeq/NR, and additional curated hydrogenase sequences.

Sequences were reviewed for subgroup consistency and completeness. A final length-based completeness filter excluded proteins shorter than 70% of their subgroup's median length.

## Repository structure

```text
MC-HYD/
├── README.md
├── LICENSE
├── LICENSE-DATA.md
├── mc_hyd_annotate.py
├── database/
│   └── MC_HYD_DB_20261025_filtered.fasta
├── hmm/
│   └── MC_HYD_DB_HMM_20261025.hmm
└── kaiju/
    ├── <MC-HYD index>.fmi
    ├── nodes.dmp
    ├── names.dmp
    └── <subgroup mapping>.tsv
```

The HMM library and Kaiju reference database are distributed separately, allowing users to install and use either annotation approach independently.

The Kaiju index and accompanying classification files must be generated from the same MC-HYD reference release.

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/MC-HYD.git
cd MC-HYD
```

### HMMER

Protein annotation requires Python 3 and HMMER 3.

On Ubuntu/Debian:

```bash
sudo apt install python3 hmmer
```

Alternatively, with Conda:

```bash
conda install -c bioconda hmmer
```

The MC-HYD Python annotation script uses only the Python standard library.

### Kaiju

Install Kaiju using Conda or Mamba:

```bash
conda install -c bioconda kaiju
```

Alternatively, install Kaiju from its official repository:

https://github.com/bioinformatics-centre/kaiju

See the official Kaiju documentation for platform-specific installation instructions.

## 1. Protein and MAG annotation using HMMER

MC-HYD provides 39 subgroup-specific profile HMMs constructed from the complete curated reference database.

The HMMs can identify hydrogenase catalytic subunits in:

- Predicted proteins from isolate genomes
- Metagenome-assembled genomes (MAGs)
- Assembled metagenomic protein catalogs
- Individual protein sequences

### Run the annotation pipeline

Input must be an amino-acid FASTA file.

```bash
python mc_hyd_annotate.py proteins.faa -o results/my_MAG
```

The script searches all input proteins against the MC-HYD HMM library using `hmmscan`.

For each protein, it:

1. Evaluates matches against all 39 subgroup-specific HMMs.
2. Selects the HMM with the highest bit score.
3. Retains the prediction if the score is at least **400 bits**.
4. Assigns the corresponding hydrogenase family and subgroup.
5. Generates a complete list of accepted hydrogenases and a subgroup summary.

Selenocysteine (`U`) is normalized to cysteine (`C`), consistent with the HMM construction procedure.

### HMM detection threshold

The validated operational detection threshold is:

**Bit score ≥ 400**

This threshold was evaluated against internal hydrogenase references, independent hydrogenase sequences, challenging non-hydrogenase homologs, the mi-faser GS-24-all curated enzyme collection, and UniProtKB/Swiss-Prot.

The threshold is intended primarily for complete hydrogenase catalytic subunits. Short or partial proteins may score below the detection threshold.

### Output files

The annotation script generates three files.

**`hydrogenases.tsv`**

Complete list of detected hydrogenases:

| protein | group | subgroup | bitscore | evalue |
|---|---|---|---:|---|
| protein_001 | [NiFe] | 1d | 874.2 | 1.2e-260 |
| protein_002 | [FeFe] | A1 | 645.3 | 4.5e-190 |
| protein_003 | [NiFe] | 4e | 589.4 | 2.1e-177 |

*Example output; values are illustrative.*

**`hydrogenase_summary.tsv`**

Number of detected hydrogenases by group and subgroup:

| group | subgroup | count |
|---|---|---:|
| [NiFe] | 1d | 1 |
| [NiFe] | 4e | 1 |
| [FeFe] | A1 | 1 |

**`hmmscan.tbl`**

Original HMMER results, including secondary matches and alternative subgroup assignments.

This file is retained to permit further inspection and alternative filtering strategies.

### Direct HMMER usage

The HMM library can also be used without the Python script.

First, index the library:

```bash
hmmpress hmm/MC_HYD_DB_HMM_20261025.hmm
```

Search an amino-acid FASTA:

```bash
hmmscan \
  --cpu 8 \
  --noali \
  --tblout results.tbl \
  -E 1e-5 \
  --domE 1e-5 \
  hmm/MC_HYD_DB_HMM_20261025.hmm \
  proteins.faa
```

The HMMER table contains raw scores and alternative matches. The 400-bit cutoff and highest-scoring-subgroup classification must be applied separately when using `hmmscan` directly.

For consistent preprocessing and output formatting, use `mc_hyd_annotate.py`.

## 2. Metagenomic read annotation using Kaiju

MC-HYD also provides a Kaiju-compatible database for direct hydrogenase classification from metagenomic sequencing reads.

Kaiju performs protein-level sequence matching without requiring metagenome assembly or ORF prediction.

The MC-HYD Kaiju database uses functional hydrogenase subgroup identifiers instead of conventional organism taxonomy.

### Database files

The `kaiju/` directory contains:

- The MC-HYD FM-index (`.fmi`)
- `nodes.dmp`
- `names.dmp`
- Subgroup identifier mapping

The index and accompanying taxonomy files must be used together.

In the examples below, `<MC-HYD index>.fmi` represents the actual FM-index filename included in the release.

### Paired-end reads

Run Kaiju on paired-end metagenomic reads:

```bash
kaiju \
  -t kaiju/nodes.dmp \
  -f "kaiju/<MC-HYD index>.fmi" \
  -i sample_R1.fastq.gz \
  -j sample_R2.fastq.gz \
  -z 8 \
  -v \
  -o sample.kaiju.tsv
```

### Single-end reads

```bash
kaiju \
  -t kaiju/nodes.dmp \
  -f "kaiju/<MC-HYD index>.fmi" \
  -i sample.fastq.gz \
  -z 8 \
  -v \
  -o sample.kaiju.tsv
```

The `-z` parameter controls the number of CPU threads.

The `-v` option provides additional match information, including alignment scores, for subsequent filtering and analysis.

### Convert Kaiju assignments to subgroup names

Kaiju uses numeric classification identifiers. The accompanying MC-HYD mapping associates these identifiers with hydrogenase subgroups.

Subgroup names can be appended using Kaiju's built-in utility:

```bash
kaiju-addTaxonNames \
  -t kaiju/nodes.dmp \
  -n kaiju/names.dmp \
  -i sample.kaiju.tsv \
  -o sample.named.tsv
```

The resulting assignments can be aggregated to calculate:

- Hydrogenase-associated read counts
- Subgroup presence and absence
- Relative subgroup distributions
- Differences in hydrogenase functional profiles across metagenomes

Because MC-HYD uses functional subgroup identifiers, these classifications should be interpreted as hydrogenase functional annotations rather than organism taxonomic assignments.

Read assignments are not direct measurements of hydrogenase activity, protein expression, or microbial cell abundance.

## 3. Validation and performance

MC-HYD HMMs and the Kaiju read-based approach were evaluated against curated positive and negative controls.

### HMM validation

The 39 subgroup HMMs were constructed from all 7,134 proteins in the filtered reference database.

The validation included:

- Internal reference recognition and cross-subgroup classification
- Independent hydrogenase sequences
- Curated non-hydrogenase homologs
- mi-faser GS-24-all
- UniProtKB/Swiss-Prot

#### HMM performance summary

| Benchmark | Result |
|---|---|
| Internal hydrogenase detection, before score filtering | 7,134/7,134 (100%) |
| Internal exact subgroup classification | 98.44% |
| Internal sensitivity at 400 bits | 99.86% |
| External hydrogenase detection, before score filtering | 372/372 (100%) |
| External exact subgroup classification | 370/372 (99.46%) |
| External sensitivity at 400 bits | 99.73% |
| Confounder specificity at 400 bits | 422/422 (100%) |
| GS-24-all catalytic-subunit detection | 9/9 (100%) |
| GS-24-all non-target rejection | 12,482/12,482 (100%) |

The two external classification discrepancies involve sequences historically annotated as [Fe] hydrogenases but assigned by the classifiers to [FeFe] group A. These remain under investigation as potential historical annotation inconsistencies.

Internal classification discrepancies were predominantly concentrated within [FeFe] group A, consistent with known difficulties distinguishing closely related group A subtypes.

### UniProtKB/Swiss-Prot screening

The MC-HYD HMM library was screened against **576,073 reviewed UniProtKB/Swiss-Prot proteins** using the 400-bit threshold.

Only **59 proteins (0.010%)** exceeded the threshold.

Of these:

- 57 had hydrogenase-related annotations.
- Two were annotated as uncharacterized proteins already represented in MC-HYD.
- 34 were exact sequence matches to the MC-HYD reference database.

The annotation review identified no unrelated-protein false positives and no rejected complete catalytic hydrogenases. Three rejected proteins described as hydrogenase catalytic subunits were explicitly annotated as incomplete sequence fragments.

These results support high specificity across a broad protein background. The Swiss-Prot annotation screening does not constitute an exhaustive estimate of sensitivity.

### Read-based validation

Kaiju was evaluated alongside three alternative functional annotation methods:

- mi-faser
- kallisto
- FunProfiler

The benchmark included internal simulated reads, independent hydrogenases, non-hydrogenase confounders, and environmental metagenomes.

Kaiju was selected for the MC-HYD read-based workflow based on its combination of:

- High read-level classification accuracy
- Low subgroup assignment ambiguity
- Specificity against challenging homologs
- Consistency with other methods across environmental samples
- Computational efficiency

The same MC-HYD subgroup classification is used across the protein-based and read-based annotation approaches.

## 4. Database provenance and reproducibility

MC-HYD incorporates publicly available hydrogenase catalytic-subunit sequences from several sources, including:

- The original HydDB database
- NCBI RefSeq and NR
- Additional published hydrogenase reference sequences

The database was expanded and curated to improve representation of hydrogenase functional diversity, including recently described and underrepresented subgroups.

The final reference database contains 7,134 sequences following completeness filtering.

The supplied profile HMMs were constructed using MAFFT for multiple sequence alignment and HMMER for profile construction.

For reproducibility, analyses should use matching versions of the reference FASTA, HMM library, and Kaiju index.

## 5. Citation

If you use MC-HYD, please cite the associated database and validation manuscript once published.

**MC-HYD database and annotation toolkit:** Manuscript in preparation.

Please also cite the original HydDB publication:

Søndergaard, D., Pedersen, C. N. S., & Greening, C. (2016). HydDB: A web tool for hydrogenase classification and analysis. *Scientific Reports*, 6, 34212. https://doi.org/10.1038/srep34212

For HMM-based annotation:

Eddy, S. R. (2011). Accelerated profile HMM searches. *PLoS Computational Biology*, 7, e1002195. https://doi.org/10.1371/journal.pcbi.1002195

For read-based annotation:

Menzel, P., Ng, K. L., & Krogh, A. (2016). Fast and sensitive taxonomic classification for metagenomics with Kaiju. *Nature Communications*, 7, 11257. https://doi.org/10.1038/ncomms11257

## 6. License and data attribution

MC-HYD is distributed as an open scientific resource.

### Software

The MC-HYD Python annotation scripts are released under the **MIT License**. See `LICENSE`.

### Database and models

Original MC-HYD curation, subgroup annotations, database organization, and derived HMM profiles are released under the **Creative Commons Attribution 4.0 International (CC BY 4.0) License**, to the extent that these contributions are subject to copyright or database rights.

See `LICENSE-DATA.md`.

### Source sequence attribution

MC-HYD incorporates publicly available protein sequences from HydDB, NCBI RefSeq/NR, and other published sources.

The original HydDB publication is available under CC BY 4.0. NCBI generally permits reuse and redistribution of its molecular sequence data, subject to any applicable third-party rights.

Source sequences retain any rights and conditions associated with their original distribution. The MC-HYD license does not supersede such rights.

Users redistributing or modifying MC-HYD should preserve database provenance and appropriately acknowledge both MC-HYD and the original sequence sources.

### Third-party software

Kaiju and HMMER are independent software packages distributed under their respective licenses. Their licenses are not modified by the MC-HYD license.

## 7. Contact

**Donato Giovannelli**  
Department of Biology  
University of Naples Federico II, Italy

Website: https://www.donatogiovannelli.com

For questions, issues, or suggestions, please use the GitHub Issues section of this repository.
