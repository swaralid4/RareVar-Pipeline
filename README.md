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

