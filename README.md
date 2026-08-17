# Local Ancestry Inference in U.S. Honey Bees

Repository containing the manuscript, data, analysis code, figures, tables, and supporting materials for:

**Fikere et al. (2026). Local ancestry inference and genomic patterns of ancestry in U.S. honey bee populations.**

Manuscript currently under **second-round review** at *Ecology and Evolution*.

## Repository Organization

The repository is organized by major analysis and reproducibility components:

```text
├── RFMix/
│   └── Local ancestry inference and ancestry tract analyses
│
├── TreeMix/
│   └── Population structure, covariance, and migration-edge analyses
│
├── code_figures/
│   └── Scripts used to generate manuscript figures
│
├── code_tables/
│   └── Scripts used to generate manuscript and supplementary tables
│
└── DataQC/
    └── Genotype quality control, filtering, and data preparation
```

### `RFMix/`

Scripts and results related to local ancestry inference, ancestry proportions, ancestry tracts, tract lengths, and ancestry-associated genomic regions.

### `TreeMix/`

TreeMix analyses of population relationships, covariance patterns, and inferred migration edges.

### `code_figures/`

Code used to generate the main and supplementary figures presented in the manuscript.

### `code_tables/`

Code used to generate the main and supplementary tables.

### `DataQC/`

Genotype quality-control, filtering, and data-preparation procedures performed before downstream analyses.

## Analytical Workflow

The main workflow consisted of:

1. Genotype quality control and variant filtering
2. Preparation of reference and target genotype datasets
3. Local ancestry inference using RFMix
4. Estimation of genome-wide ancestry proportions and ancestry tracts
5. Analysis of geographic patterns of ancestry and tract length
6. Identification of ancestry-associated genomic regions
7. Gene identification and Gene Ontology enrichment
8. TreeMix analysis of population relationships and covariance
9. Generation of manuscript figures and tables

## Key Methods

* **Local ancestry:** RFMix v2
* **Population structure:** TreeMix
* **Genomic region analysis:** BEDTools
* **GO enrichment:** `topGO`
* **Statistical analysis and visualization:** R

Specific software versions and analysis parameters are described in the manuscript Methods.

## Manuscript Status

**Journal:** *Ecology and Evolution*
**Status:** Second-round review

The repository contains the analysis materials and supporting files associated with the current manuscript revision.
