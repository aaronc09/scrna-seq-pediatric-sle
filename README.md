# LIANA+ Inferred Cell-Cell Interaction Analysis Reveals Differences in SLE-Associated Immune Cell Communication Between Pediatric and Adult Cohorts

## Overview

This repository contains the code and supporting materials for the study:

**“LIANA+ Inferred Cell-Cell Interaction Analysis Reveals Differences in SLE-Associated Immune Cell Communication Between Pediatric and Adult Cohorts.”**

**Authors:** Aaron Choi and Courtney Hatton

This study reanalyzed publicly available single-cell RNA-sequencing (scRNA-seq) data to investigate differences in systemic lupus erythematosus (SLE)-associated immune-cell communication between pediatric and adult cohorts.

## Purpose

The purpose of this study was to investigate whether SLE-associated changes in immune-cell communication differ between pediatric and adult patients. Using publicly available single-cell RNA-sequencing data, we compared inferred cell-cell communication patterns between SLE and healthy donors in pediatric and adult cohorts.

The primary comparison was:

**(Pediatric SLE − Pediatric Healthy) − (Adult SLE − Adult Healthy)**

This approach evaluates whether the SLE-versus-healthy difference in inferred communication differs between pediatric and adult cohorts rather than directly comparing pediatric SLE with adult SLE.

## Outcome

The analysis identified **220 ligand-receptor and cell-type combinations** in which the SLE-versus-healthy difference in inferred communication differed significantly between the pediatric and adult cohorts at FDR < 0.05.

Non-classical monocytes were the most frequently represented receiving cell type among the significant results, appearing in **85 of 530 eligible combinations (16.0%)**.

Sensitivity analyses showed that some findings depended on analytical choices and technical batch. Therefore, the results are best viewed as exploratory findings that identify candidate immune-cell communication patterns for future validation rather than established age-specific biological mechanisms.

## Dataset

The analysis used publicly available peripheral blood mononuclear cell (PBMC) single-cell RNA-sequencing data from:

**NCBI Gene Expression Omnibus (GEO): GSE135779**

The final analysis included:

- **56 donors**
- 33 pediatric SLE donors
- 11 pediatric healthy donors
- 7 adult SLE donors
- 5 adult healthy donors
- **298,610 retained cells**

No new patient data were collected for this study.

## Analysis

The workflow included:

- Single-cell quality control and preprocessing with Scanpy
- Doublet detection with Scrublet
- Donor-aware integration with Harmony
- Leiden clustering and immune-cell annotation
- CellTypist-supported annotation review
- Donor-level cell-cell communication inference using LIANA+
- Ligand-receptor interaction analysis using a consensus resource
- Disease-by-cohort interaction modeling
- HC3 heteroskedasticity-robust standard errors
- Benjamini-Hochberg false discovery rate correction
- Multiple sensitivity and robustness analyses

The **donor**, rather than the individual cell, was treated as the independent biological replicate.

## Main Findings

Among **2,335 estimable ligand-receptor and sending-receiving cell-type combinations**, 220 met FDR < 0.05 in the primary cohort-by-SLE analysis.

Of these:

- 3 interaction coefficients were positive
- 217 interaction coefficients were negative
- Non-classical monocytes were the receiving cell type in 85 of 530 eligible combinations (16.0%)

Sensitivity and robustness analyses evaluated:

- Female-only donors
- Alternative cell-annotation choices
- Alternative ligand-receptor resources
- Communication-score availability
- Low-expression interactions
- Donor influence
- Technical batch effects

The findings showed sensitivity to several analytical choices. Technical batch was also substantially imbalanced across the pediatric and adult cohorts. These limitations prevent attributing the observed differences specifically to biological age.

## Repository Contents

This repository contains code and supporting materials for:

1. Data preprocessing and quality control
2. Cell clustering and annotation
3. LIANA+ cell-cell communication inference
4. Donor-level statistical modeling
5. Sensitivity and robustness analyses
6. Figure generation
7. Supplemental analyses and tables

Complete primary and sensitivity-analysis results are available in:

`Results/supplemental_tables/`

## Reproducibility

The primary analysis was performed in **Python 3.12.10**.

Major packages included:

- Scanpy 1.12.3
- AnnData 0.13.2
- LIANA 1.8.1
- CellTypist 1.7.1
- Scrublet 0.2.3
- harmonypy 0.0.10
- leidenalg 0.12.0
- igraph 1.0.0
- umap-learn 0.5.12
- NumPy 2.5.2
- pandas 2.3.3
- SciPy 1.18.0
- statsmodels 0.14.6

LIANA 1.8.1 was used with the documented process-local correction to the `magnitude_rank` score-column ranking described in the manuscript and analysis code.

## Data Availability

The original single-cell RNA-sequencing dataset is publicly available through the **NCBI Gene Expression Omnibus (GEO)**:

**GSE135779**

Supplemental tables containing the complete primary and sensitivity-analysis results are available within this repository.

## Interpretation and Limitations

LIANA+ infers potential cell-cell communication from RNA expression and ligand-receptor resources. These predictions do not directly measure:

- Protein abundance
- Ligand-receptor binding
- Physical cell-cell contact
- Pathway activation
- Absolute signaling intensity

This study uses a single cross-sectional blood dataset and does not include an independent replication cohort.

Technical batch was also imbalanced between pediatric and adult donors. Although sensitivity analyses were performed, the available study design cannot fully separate biological age from technical and cohort-related effects.

Therefore, the results should be interpreted as **exploratory cohort-associated differences** rather than evidence of causal mechanisms, confirmed age-specific signaling changes, or therapeutic targets.

Independent, technically balanced cohorts and protein-level and functional validation are needed to evaluate these findings further.

## Citation

If you use the code or results from this repository, please cite the associated manuscript:

**Choi A, Hatton C. LIANA+ Inferred Cell-Cell Interaction Analysis Reveals Differences in SLE-Associated Immune Cell Communication Between Pediatric and Adult Cohorts.**

Publication information and DOI will be added when available.

## Authors

**Aaron Choi**  
Northern Valley Regional High School at Old Tappan

**Courtney Hatton**  
UMass Chan School of Medicine
