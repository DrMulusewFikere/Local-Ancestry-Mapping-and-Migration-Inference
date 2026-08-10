# Local Ancestry Inference in U.S. Honey Bees

## Project overview

This repository contains the manuscript, supplementary materials, analysis scripts,
and supporting results for the study:

**Fikere et al. (2026). Local ancestry inference and genomic patterns of ancestry
in U.S. honey bee populations.**

The manuscript is currently under **second-round review** at *Ecology and Evolution*.

This repository contains the final materials prepared for the current revision.
Files corresponding to the final submission are identified with the suffix:

`_Final`

---

## Repository organization

### `manuscript/`

Contains the manuscript and submission-related documents.

- Main revised manuscript
- Response to reviewers
- Cover letter

Files with the `_Final` suffix correspond to the versions prepared for submission.

---

### `figures/`

Contains the figures used in the manuscript.

#### `figures/main/`

Main-text figures.

#### `figures/supplementary/`

Supplementary figures, including:

- Local ancestry analyses
- Regression analyses
- Gene Ontology enrichment
- TreeMix analyses
- Additional population-genetic analyses

---

### `tables/`

Contains tables associated with the manuscript.

#### `tables/main/`

Tables included in the main manuscript.

#### `tables/supplementary/`

Supplementary tables, including:

- TreeMix covariance results
- Gene Ontology enrichment results for A-, C-, M-, and O-lineage ancestry

---

### `data/`

Contains input data and reference files used in the analyses.

#### `data/metadata/`

Sample metadata, including geographic and sampling information.

#### `data/genotype/`

Genotype and variant datasets used for downstream analyses.

#### `data/local_ancestry/`

Input and output data associated with RFMix local ancestry inference.

#### `data/reference/`

Reference genome and annotation files used in the analyses.

---

### `analysis/`

Contains analysis outputs and intermediate results organized by analytical component.

#### `analysis/local_ancestry/`

Local ancestry inference and ancestry proportion analyses.

#### `analysis/regression/`

Regression analyses examining relationships between ancestry proportions,
latitude, longitude, and other geographic variables.

#### `analysis/haplotype_tracts/`

RFMix ancestry tract identification and tract-length analyses.

#### `analysis/TreeMix/`

TreeMix analyses of population relationships, covariance, and inferred migration
edges.

#### `analysis/GO_enrichment/`

Gene Ontology enrichment analyses for genomic regions associated with A-, C-, M-,
and O-lineage ancestry.

---

### `scripts/`

Contains analysis scripts used to generate the results.

Scripts are organized according to the corresponding analytical component:

- `local_ancestry/`
- `regression/`
- `haplotype_tracts/`
- `TreeMix/`
- `GO_enrichment/`

---

### `results/`

Contains analysis results and intermediate files generated during the study.

Results are organized according to the corresponding analysis:

- `local_ancestry/`
- `regression/`
- `haplotype_tracts/`
- `TreeMix/`
- `GO_enrichment/`

---

## Major analytical workflow

The overall analysis consisted of the following major steps:

1. Genomic data generation and quality control
2. Genotype imputation and variant filtering
3. Preparation of reference and target genotype datasets
4. Local ancestry inference using RFMix
5. Estimation of genome-wide ancestry proportions
6. Analysis of geographic patterns of ancestry
7. Identification and characterization of ancestry tracts
8. Analysis of ancestry tract length and geographic variation
9. TreeMix analysis of population relationships and covariance
10. Identification of ancestry-associated genomic regions
11. Gene identification within ancestry-associated regions
12. Gene Ontology enrichment analysis
13. Integration and visualization of ancestry-specific genomic patterns

---

## Ancestry-associated genomic regions

Regions with elevated ancestry probabilities were identified using the specified
regional SNP/ancestry criterion described in the manuscript.

Candidate genomic regions were intersected with the *Apis mellifera*
Amel_HAv3.1 genome annotation using BEDTools to identify overlapping genes.

These candidate genes were subsequently used for Gene Ontology enrichment analysis.

---

## Gene Ontology enrichment

GO enrichment was performed using `topGO` with:

- Biological Process (BP) ontology
- `weight01` algorithm
- Fisher's exact test
- `nodeSize = 10`
- Benjamini-Hochberg false-discovery-rate correction

Separate analyses were performed for A-, C-, M-, and O-lineage-associated
candidate genes.

---

## TreeMix analysis

TreeMix was used to evaluate population relationships, covariance patterns, and
inferred migration edges among geographic population groups.

The final TreeMix analysis used five geographic clusters and four migration edges.

Migration edges are interpreted as covariance patterns consistent with gene flow
rather than as definitive evidence of contemporary migration between individual
states. Shared ancestry and historical demographic processes may also contribute
to the observed patterns.

The TreeMix results are presented in the supplementary figures and tables.

---

## Haplotype tract analysis

RFMix-inferred ancestry tracts were used to evaluate the geographic distribution
and physical length of ancestry-specific haplotype tracts.

Physical tract length was calculated from the genomic coordinates of each tract:

`(EndPos - StartPos) / 1,000,000`

resulting in tract length measured in megabases (Mb).

The relationship between A-lineage tract length and latitude was evaluated using
linear regression.

---

## Software

Major software and R packages used in the analyses include:

- R
- RFMix v2
- TreeMix
- BEDTools
- topGO
- ggplot2
- dplyr
- tidyr
- GenomicRanges
- patchwork
- ggpmisc

Specific software versions and analysis parameters are reported in the manuscript
Methods section where applicable.

---

## Manuscript status

Journal: *Ecology and Evolution*

Manuscript status: **Second-round review**

The repository corresponds to the materials prepared for the current revision.
