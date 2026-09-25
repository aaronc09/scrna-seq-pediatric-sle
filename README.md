# Pediatric vs. Adult SLE Cell-Cell Communication

## LIANA+ Analysis of Single-Cell RNA-Seq Data

This repository contains the code and analysis for the study:

**“LIANA+ Inferred Cell-Cell Interaction Analysis Reveals Differences in SLE-Associated Immune Cell Communication Between Pediatric and Adult Cohorts”**

**Authors:** Aaron Choi and Courtney Hatton

---

## Overview

Systemic lupus erythematosus (SLE) is an autoimmune disease in which the immune system attacks the body's own tissues. SLE can occur in both children and adults, but pediatric and adult disease can differ in clinical and immune characteristics.

Most single-cell studies examine which genes are expressed by different immune-cell populations. In this project, we asked a related question:

> **Do SLE-associated patterns of communication between immune cells differ between pediatric and adult cohorts?**

To investigate this, we reanalyzed publicly available single-cell RNA-sequencing (scRNA-seq) data and used **LIANA+** to infer potential cell-cell communication based on the expression of known ligand-receptor pairs.

Importantly, these analyses identify **candidate communication patterns from RNA expression**. They do not directly demonstrate physical cell-cell interactions, protein binding, or signaling activity.

---

## Study Design

The analysis used publicly available peripheral blood mononuclear cell (PBMC) scRNA-seq data from:

**GEO accession: GSE135779**

The analyzed dataset included **56 donors**:

| Group | Donors |
|---|---:|
| Pediatric SLE | 33 |
| Pediatric healthy | 11 |
| Adult SLE | 7 |
| Adult healthy | 5 |
| **Total** | **56** |

After quality control and revised cell-type annotation, the primary analysis included:

- **233,770 pediatric cells**
- **64,840 adult cells**
- **298,610 cells total**

Major immune-cell populations included:

- CD4 T cells
- CD8 T cells
- B cells
- NK cells
- Classical monocytes
- Non-classical monocytes
- Plasma cells
- Plasmacytoid dendritic cells (pDCs)

---

## Main Research Question

Simply comparing pediatric SLE with adult SLE could be misleading because healthy children and healthy adults may already have biological differences.

Instead, we first compared SLE with healthy donors **within each age cohort**, and then compared those changes:

```text
(pediatric SLE - pediatric healthy)
              -
(adult SLE - adult healthy)
```

This is the **age-group-by-SLE interaction** used in the primary statistical analysis.

In simpler terms, the analysis asked:

> **Is the change in inferred immune-cell communication associated with SLE different in the pediatric cohort than in the adult cohort?**

---

## Analysis Workflow

The overall workflow was:

```text
GSE135779 scRNA-seq data
        |
        v
Quality control
        |
        v
Doublet detection with Scrublet
        |
        v
Normalization and log transformation
        |
        v
Highly variable gene selection
        |
        v
PCA
        |
        v
Harmony integration by donor
        |
        v
Leiden clustering + UMAP
        |
        v
Cell-type annotation
        |
        v
Donor-level LIANA+ analysis
        |
        v
Ligand-receptor communication scores
        |
        v
Age-group x SLE statistical model
        |
        v
Multiple-testing correction
        |
        v
Sensitivity analyses
```

A major feature of the study is that **each donor, rather than each individual cell, was treated as the independent biological replicate**.

This avoids treating thousands of cells obtained from the same person as if they were thousands of independent patients.

---

## Cell-Cell Communication Analysis

Cell-cell communication was inferred using **LIANA+ consensus rank aggregation**.

LIANA+ uses known ligand-receptor relationships together with single-cell gene-expression data to identify potential communication between cell populations.

Conceptually:

```text
Source cell
    |
    v
Ligand
    |
    v
Receptor
    |
    v
Target cell
```

LIANA+ inference was performed **separately for each donor**.

A cell type was included for a donor only when at least **10 cells** of that type were available.

Ligand-receptor interactions were retained only when the ligand and receptor were detected in at least **10% of cells** in their respective source and target populations.

The LIANA+ magnitude rank was converted to:

```text
communication score = 1 - magnitude_rank
```

so that larger values represented stronger inferred communication.

---

## Statistical Analysis

For each eligible ligand-receptor and source-target cell combination, the following donor-level model was fitted:

```text
communication score =
    disease status
    + age group
    + disease status x age group
```

The interaction coefficient represented:

```text
(pediatric SLE - pediatric healthy)
              -
(adult SLE - adult healthy)
```

Only interactions with data from at least **5 donors in each of the four groups** were tested.

The analysis used:

- HC3 robust standard errors
- 95% confidence intervals
- Benjamini-Hochberg false discovery rate (FDR) correction

The primary significance threshold was **FDR < 0.05**.

---

## Main Results

The primary analysis tested **2,374 ligand-receptor/source-target combinations**.

Of these:

- **171 interactions had FDR < 0.05**
- **239 interactions had FDR < 0.10**
- **37** significant interactions had positive interaction coefficients
- **134** had negative interaction coefficients

These results indicate that SLE-associated changes in inferred immune-cell communication were not identical between the pediatric and adult cohorts.

---

## Non-Classical Monocytes Were Prominent

One of the clearest overall patterns involved **non-classical monocytes**.

Among the 171 significant interactions, non-classical monocytes were the receiving cell type in **57 interactions**, including:

- **12** from CD4 T cells
- **9** from CD8 T cells

Classical monocytes were also frequent communication sources, contributing **27 significant interactions**, including 9 directed toward non-classical monocytes.

These counts do **not** prove that non-classical monocytes are more active, more important, or causal in SLE. They show that monocyte-associated communication was prominent among the statistically significant differences identified in this analysis.

---

## Some Communication Patterns Changed in Opposite Directions

The pediatric-adult differences were not caused by a simple overall increase or decrease in inferred communication.

Several leading interactions changed in **opposite directions** between the pediatric and adult cohorts.

Examples that increased with SLE in the pediatric cohort but decreased in the adult cohort included:

- **TIMP2-ITGB1**
- **CD48-CD244**
- **SELPLG-ITGB2**

Other interactions generally showed the opposite pattern, including:

- **TGFB1-ITGB1**
- **CD48-CD2**
- **CD52-SIGLEC10**

This suggests that pediatric and adult SLE may differ in the **pattern and direction** of specific inferred immune-cell communication changes rather than simply in the overall amount of communication.

---

## Sensitivity Analyses

Several sensitivity analyses were performed to examine how robust the findings were to analytical choices.

### Female-Only Analysis

Because all analyzed adult donors were female while the pediatric cohort included four male donors, the primary analysis was repeated using only female donors.

Results were highly similar to the primary analysis:

**Pearson r = 0.995**

The female-only analysis identified:

- **145 interactions at FDR < 0.05**
- **220 interactions at FDR < 0.10**
- **141 of the 171** primary significant interactions remained significant

This suggests that the overall coefficient pattern was not primarily driven by the sex imbalance between cohorts.

### Cell-Type Annotation Sensitivity

The revised primary annotation was compared with the original annotation that excluded the full T-rich clusters.

The original-exclusion analysis identified:

- **123 significant interactions**
- **2,527 interactions tested**

Across 2,371 tests shared between the two analyses:

**Pearson r = 0.910**

### Ligand-Receptor Resource Sensitivity

The primary LIANA+ consensus resource was also compared with the CellPhoneDB resource bundled with LIANA.

Using the revised annotations with CellPhoneDB:

- **514 interactions were tested**
- **29 met FDR < 0.05**

For shared tests, coefficients correlated with the primary consensus analysis at:

**Pearson r = 0.828**

Importantly, **non-classical monocytes remained the most frequent significant receiving cell type** across the primary analysis and these sensitivity analyses.

Individual ligand-receptor findings, however, depended more strongly on annotation and resource choices.

---

## How to Interpret the Findings

This study uses computational inference.

LIANA+ combines RNA-expression measurements with databases of known or proposed ligand-receptor relationships.

Therefore, a significant interaction in this study means:

> **The RNA-based inferred communication pattern differed in its SLE-associated change between the pediatric and adult cohorts.**

It does **not** by itself demonstrate:

- Physical communication between the cells
- Ligand-receptor protein binding
- Protein abundance
- Downstream pathway activation
- Causation
- A validated drug target

Protein-level studies, functional experiments, and independent cohorts would be needed to confirm the biological importance of individual interactions.

---

## Important Limitations

### Sample Size

Although almost 300,000 cells were analyzed, the biological sample consisted of **56 donors**, including only 12 adults.

### Age and Technical Batch

The pediatric and adult cohorts differed in technical batch, so biological age cannot be completely separated from other cohort differences.

### Clinical Confounding

Disease activity, medications, disease duration, ancestry, and clinical manifestations could not all be controlled.

### Sex Imbalance

All analyzed adults were female, although the female-only sensitivity analysis produced highly similar results.

### Computational Inference

LIANA+ predicts candidate communication patterns from RNA expression rather than directly measuring protein-level signaling.

### Annotation and Database Uncertainty

Individual results can depend on how cells are annotated and which ligand-receptor resource is used.

### Peripheral Blood Only

The analysis used PBMCs. Communication patterns within affected tissues such as the kidney or skin may differ from those observed in blood.

---

## Why This Study Matters

Immune cells do not function independently. They communicate through signaling molecules, receptors, direct contact, and complex signaling networks.

Understanding these communication networks may provide information that cannot be obtained by studying individual genes or cell types alone.

This may be especially relevant as more treatments are designed to target specific immune cells, receptors, cytokines, and signaling pathways.

The findings from this study suggest that some SLE-associated immune communication patterns may differ between pediatric and adult cohorts.

However, the results should be considered **hypothesis-generating**. Further studies are needed to determine whether these communication differences are reproducible, biologically functional, related to disease severity or organ involvement, or relevant to treatment response.

---

## Data Availability

The scRNA-seq data used in this study are publicly available through the NCBI Gene Expression Omnibus:

**GSE135779**

The study is a secondary analysis of deidentified, publicly available data.

---

## Software and Major Tools

The analysis used tools including:

- Python 3.12
- Scanpy
- AnnData
- Scrublet
- Harmony
- Leiden clustering
- UMAP
- CellTypist
- LIANA+
- Python statistical analysis tools

See the analysis scripts and environment files in this repository for the exact computational workflow and package versions.

---

## Reproducibility

The repository contains the computational workflow used for:

- Data preprocessing
- Quality control
- Doublet detection
- Dimensionality reduction
- Cell-type annotation
- Annotation review
- LIANA+ cell-cell communication inference
- Donor-level statistical analysis
- Sensitivity analyses
- Table and figure generation

Intermediate and final outputs are organized so that the major steps of the analysis can be inspected and reproduced.

---

## Citation

**Choi A, Hatton C.**  
*LIANA+ Inferred Cell-Cell Interaction Analysis Reveals Differences in SLE-Associated Immune Cell Communication Between Pediatric and Adult Cohorts.*

---

## Authors

**Aaron Choi**  
Northern Valley Regional High School at Old Tappan

**Courtney Hatton**  
UMass Chan School of Medicine
