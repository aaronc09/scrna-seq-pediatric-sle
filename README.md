# LIANA+ Inferred Cell-Cell Interaction Analysis Reveals Differences in SLE-Associated Immune Cell Communication Between Pediatric and Adult Cohorts

**Authors:** Aaron Choi and Courtney Hatton

This repository contains analysis code and results for a secondary analysis of GSE135779 peripheral-blood single-cell RNA-sequencing data. The analysis includes 56 donors (33 pediatric SLE, 11 pediatric healthy, 7 adult SLE, and 5 adult healthy) and 298,610 retained cells.

## Current results and supplemental tables

- [Supplemental Tables S1-S9](Results/supplemental_tables/): numbered publication files, with panel descriptions in the README.
- [Figure and table guide](Results/reviewer_revision/figure_table_numbering.txt): Figures 1-5, main Tables 1-2, and Supplemental Tables S1-S9.
- [Current analysis guide](Results/reviewer_revision/README.txt): corrected inference, statistical models, verification, and sensitivity analyses.
- [All primary models](Results/reviewer_revision/primary_models.csv) and [220 significant combinations](Results/reviewer_revision/supplementary_all_significant.csv).
- [Publication figures](Results/manuscript_figures/) and [main table assets](Results/manuscript_tables/).
- [Editable Table 2](Results/reviewer_revision/table2_editable.html), with four illustrative significant combinations.

The primary analysis found 220 significant combinations among 2,335 estimable tests at FDR < 0.05. Three interaction coefficients were positive and 217 were negative. The contrast is:

```text
(pediatric SLE - pediatric healthy) - (adult SLE - adult healthy)
```

Each donor is an independent biological replicate. The outcome is a relative LIANA magnitude-rank score, analyzed with HC3 robust errors and Benjamini-Hochberg correction. Results are exploratory cohort associations: batch imbalance, annotation, resource choice, and score availability limit biological interpretation. These scores do not measure absolute signaling or establish biological age effects.

## Reproduction

Use Python 3.12.10 and the recorded dependencies in [requirements.txt](requirements.txt) or [environment.yml](environment.yml). Obtain the public expression data from GEO accession **GSE135779**; large expression matrices and local caches are excluded from Git. The two root spreadsheets supply sample mappings and source clinical metadata used by the pipeline.

See [Notebooks/README.md](Notebooks/README.md) for stages and the [current analysis guide](Results/reviewer_revision/README.txt) for prerequisites. Earlier stages establish expression objects, annotations, per-method donor scores, and comparison outputs. The corrected workflow ends with:

```powershell
python Notebooks/18_REVIEWER_REVISION.py cached
python Notebooks/18_REVIEWER_REVISION.py infer --archive-root "C:/path/to/GSE135779_transfer"
python Notebooks/19_SUMMARIZE_REVIEWER_REVISION.py
```

LIANA 1.8.1 requires the process-local ranking correction implemented in stage 18. Running stage 17 alone produces historical, uncorrected rank-score results. Follow the prerequisite workflow before these commands; compact repository outputs do not replace the expression and donor-score inputs.

## Repository contents

| Location | Purpose |
| --- | --- |
| `Notebooks/` | Preprocessing notebooks, analysis scripts, and shared modules |
| `Results/reviewer_revision/` | Authoritative corrected models, resource snapshots, diagnostics, and verification |
| `Results/supplemental_tables/` | Numbered publication copies of supplemental data |
| `Results/manuscript_figures/`, `Results/manuscript_tables/` | Current publication assets |
| `Results/qc/`, `Results/annotation_review/`, `Results/provenance/` | QC, annotation evidence, and provenance |
| Other `Results/` directories | Earlier workflow inputs and historical comparison results |

Historical outputs remain where scripts expect them because they support correction comparisons and reproducibility. Their rank-score results are superseded by `Results/reviewer_revision/`. Numbered supplemental copies are intentionally retained alongside the source analysis files so manuscript links are stable and readable.

Local manuscript drafts, reviewer correspondence, submission exports, caches, environments, and temporary files are excluded through `.gitignore`. Upload using Git so these exclusions are respected.

## Citation and license

Choi A, Hatton C. *LIANA+ Inferred Cell-Cell Interaction Analysis Reveals Differences in SLE-Associated Immune Cell Communication Between Pediatric and Adult Cohorts.* Publication details will be added when available.

See [LICENSE](LICENSE) for repository licensing.
