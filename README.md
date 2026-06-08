# RareVar-Pipeline
*Rare disease variant prioritization from raw VCF to ranked candidates.*

[![Python](https://img.shields.io/badge/Python-3.10+-blue?logo=python)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Data: 1000 Genomes](https://img.shields.io/badge/Data-1000%20Genomes-orange)](https://www.internationalgenome.org/)
[![Build: hg38](https://img.shields.io/badge/Genome%20Build-hg38%2FGRCh38-blue)](https://genome.ucsc.edu)

## Overview

A lightweight, reproducible bioinformatics pipeline for identifying potentially pathogenic rare variants associated with genetic disease. Built and benchmarked on an HP ProBook 440 G8 (Intel Core i5-1135G7 × 8, 16 GB RAM, Ubuntu 24.04 LTS) — full pipeline completes in under 9 minutes on a standard consumer laptop.

When a patient with a suspected rare genetic disease is sequenced, the resulting VCF file contains millions of variants. The fundamental challenge is that the causative mutation — the single variant actually responsible for the disease is buried within that noise. Manual inspection is impossible at scale.
This pipeline automates the triage process — filtering millions of variants down to a ranked shortlist of high-priority candidates using allele frequency, gene identity, and clinical significance as evidence. The same class of tool is used in hospitals and research institutes worldwide, including NHS Genomics England, the Broad Institute, and Mayo Clinic.

---

## Dependencies

```
cyvcf2 · pandas · numpy · scipy · matplotlib · seaborn · biopython · tqdm · jupyter · notebook · ipykernel
```

```bash
pip install -r requirements.txt
```

---

## Pipeline Stages

```
Raw VCF (chr17)
    │
    ▼
Stage 1 — VCF Parsing          parse_vcf.py
    │   Extract CHROM, POS, REF, ALT, QUAL, AF, GENE, CLNSIG
    ▼
Stage 2 — EDA                  notebooks/01_EDA.ipynb
    │   Inspect quality, AF distribution, missing values
    ▼
Stage 3 — Filtering            filter_variants.py
    │   QUAL > 30  ·  AF < 0.01  ·  AF not null
    ▼
Stage 4 — Gene Annotation      annotate_genes.py   ← custom addition
    │   Positional lookup against hg38 refGene coordinates
    ▼
Stage 5 — Disease Gene Match   match_genes.py
    │   Flag variants in OMIM-curated rare disease genes
    ▼
Stage 6 — ClinVar Annotation   annotate_clinvar.py
    │   Join against ClinVar hg38
    ▼
Stage 7 — Prioritization       prioritize.py
    │   Score each variant across 5 evidence dimensions
    ▼
Stage 8 — Visualization        visualize.py
        3 publication-ready figures saved to figures/
```

---
## Quick Start

```bash
# 1. Clone
git clone https://github.com/swaralid4/RareVar-Pipeline
cd RareVar-Pipeline

# 2. Create environment
conda create -n rarevar python=3.10 -y
conda activate rarevar
pip install -r requirements.txt

# 3. Download data (~500 MB)
wget -P data/raw/ \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/supporting/GRCh38_positions/ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

wget -P data/raw/ \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
gunzip -k data/raw/refGene.txt.gz

# 4. Run
python3 run_pipeline.py
```

Expected runtime: **~8 minutes** on a standard laptop.

---

## Project Structure

```
RareVar-Pipeline/
├── data/
│   ├── raw/                  # VCF + refGene (git-ignored, large files)
│   ├── processed/            # Intermediate CSVs (git-ignored)
│   └── results/              # Final prioritized output (git-ignored)
├── scripts/
│   ├── parse_vcf.py          # VCF → CSV with AF recomputation
│   ├── filter_variants.py    # AF < 0.01 filter
│   ├── build_gene_list.py    # OMIM disease gene JSON
│   ├── annotate_genes.py     # hg38 positional gene assignment
│   ├── match_genes.py        # Disease gene flagging
│   ├── annotate_clinvar.py   # ClinVar hg38 join
│   ├── prioritize.py         # Evidence scoring
│   └── visualize.py          # 3 publication figures
├── notebooks/
│   └── eda.ipynb             # Exploratory data analysis
├── figures/                  # Output plots (git-ignored)
├── docs/
│   └── index.html            # Lab notebook (GitHub Pages)
├── run_pipeline.py           # Single-command pipeline runner
├── requirements.txt
├── .gitignore
└── README.md
```

## Data Sources

| Source | Purpose |
|--------|---------|
| [1000 Genomes GRCh38 realigned (2019)](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/) | chr17 VCF with real population AF |
| [UCSC hg38 refGene](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/) | Gene coordinate annotation |
| [ClinVar hg38 VCF](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) | Pathogenicity classifications |
| [OMIM](https://omim.org/) | Disease gene curation |

---

## License

MIT © 2026 — Swarali Deshpande
