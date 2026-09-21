# Pediatric and adult SLE cell-cell communication

Research code and supporting results for an observational analysis of peripheral-blood single-cell RNA sequencing from GEO series **GSE135779**. The analysis compares inferred ligand-receptor communication across pediatric healthy, pediatric SLE, adult healthy, and adult SLE donors.

**Current primary analysis:** reviewed T-cell restoration, 298,610 cells from 56 donors, and 171 significant interactions among 2,374 tests. See [the revised primary package](Results/revised_primary/README.md) and [annotation review](Results/annotation_review/README.md). The original full-cluster exclusions are retained as a sensitivity comparison. Visible manuscript figures and tables correspond to the revised primary analysis.

## Repository contents

| Location | Contents |
|---|---|
| `Notebooks/` | Preprocessing, annotation, LIANA inference, statistical models, figure generation, and supporting-data exports |
| `Results/manuscript_tables/` | Main and supplementary tables as machine-readable CSV files |
| `Results/manuscript_figures/` | Research figures |
| `Results/qc/` | Quality-control summaries, annotation evidence, excluded-cell identifiers, and donor coverage |
| `Results/revised_primary/` | Current primary models, female-only models, exact labels, and comparison results |
| `Results/annotation_review/` | Marker, QC, doublet-score, and donor consistency review supporting the revision |
| `Results/sensitivity_analysis/` | Historical sensitivity results using the original exclusions |
| `Results/severity_analysis/` | Clinical metadata and historical original-exclusion model results |
| `Results/provenance/` | Recorded software versions, run information, and source-file reconciliation |
| `GSE135779 GSM IDs.xlsx` | Study-to-GEO sample mapping |
| `suppdata.xlsx` | Source supplementary metadata |

This repository contains research code and supporting results. Manuscript drafts and journal submission packages are maintained separately. Local analysis caches under `_local_submission/` are excluded from Git; omit that directory if uploading files through the GitHub website.

Current research figures and tables are kept only in `Results/manuscript_figures/` and `Results/manuscript_tables/`. Authoritative models remain in `Results/revised_primary/`. The 22 numbered historical run records have been combined into `Results/provenance/run_manifests.json`, keyed by their original filenames. Script 17 with `--publish` now removes duplicate presentation staging files and consolidates those run records automatically. Scientific provenance and original-exclusion comparisons are retained to support the research.

## Environment and input data

The recorded analysis environment uses **Python 3.12.10**. From the repository root on Windows:

```powershell
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\Activate.ps1
```

An environment specification is also provided in `environment.yml`; recorded run versions are under `Results/provenance/`.

Raw GEO matrices, large processed AnnData objects, and intermediate donor-level LIANA outputs are not bundled. To rerun the full analysis, obtain the GSE135779 source data, use the sample-mapping spreadsheet, and follow the input paths and preprocessing steps in notebooks 01A and 01D. The pipeline uses `GSE135779_RAW/`, `child_individual_h5ad/`, and `adult_individual_h5ad/` for local data. These directories are ignored by Git.

## Analysis order

All active project code is in [`Notebooks/`](Notebooks/README.md). The `.ipynb` notebooks, `.py` analysis scripts, and shared helper modules are complementary parts of one workflow, not duplicate versions.

The original pipeline below recreates historical full-cluster-exclusion results. Complete scripts 13–17 afterward to reproduce the current primary analysis; run script 17 with `--publish` last to refresh manuscript figures and tables.

Run notebooks from `Notebooks/` with the project environment selected:

1. `01A` through `01F`: preprocessing, cell annotation, and donor splitting.
2. `02A` through `02D`: donor-level LIANA inference.
3. `03A` and `03B`: clinical metadata and donor-level score extraction.
4. `04A` through `04F`: secondary group comparisons.
5. `05`: primary age-by-disease interaction models.
6. `06A` and `06B`: disease-severity analyses.

Run the following scripts from the repository root in this order:

```powershell
python Notebooks/07_GENERATE_BATCH_CHECK_FIGURES.py
python Notebooks/08_PUBLICATION_SENSITIVITY_ANALYSES.py
python Notebooks/09_GENERATE_MANUSCRIPT_SUMMARIES.py
python Notebooks/12_GENERATE_SUBMISSION_AUDITS.py
python Notebooks/10_GENERATE_MANUSCRIPT_UMAP.py
```

Script 12 retains its original filename and generates clinical-covariate summaries, model donor coverage, and supporting exclusion records. It can run from the compact saved results without refitting models. Scripts requiring expression matrices or intermediate scores need the preceding analysis outputs. Rerunning scripts can replace their generated outputs.

Randomized steps use fixed seeds recorded in the code: 0 for preprocessing/QC, 42 for network layouts, and deterministic interaction-specific bootstrap seeds derived from 0.

## Supplementary data

Under `Results/manuscript_tables/`:

- **S1:** `supplementary_adult_sample_exclusions.csv`.
- **S2:** `supplementary_cell_type_coverage.csv`, with full donor-by-cell-type counts in `supplementary_donor_cell_type_coverage.csv`.
- **S3:** `supplementary_all_fdr05_interactions.csv`, including contributing-donor counts for each group.

CSV files preserve statistical precision. Clinical source codes, including `ND`, are preserved; `ND` is not interpreted as absence of treatment. Recorded race and ethnicity are clinical categories, not genetic ancestry. These descriptive summaries do not adjust the primary model.

## Annotation evidence and reproducibility limits

`Results/qc/` contains cluster counts, donor-by-cluster counts, excluded-cell identifiers, full/top-20 marker rankings, and canonical-marker expression summaries. The recovered pre-exclusion objects reproduce the saved donor counts and retained-cell marker summaries. Source records and reconciliation results are in `Results/provenance/archived_annotation_audit.json`; source paths are relative to an external archive root.

To repeat that check with the external archive available:

```powershell
python Notebooks/13_CHECK_ARCHIVED_ANNOTATION.py --archive-root "C:/path/to/GSE135779_transfer"
```

This optional script requires pandas, numpy, and h5py. It expects the cohort AnnData files in the `child_individual_h5ad/` and `adult_individual_h5ad/` subdirectories and marker tables in their original notebook subdirectories under `Notebook_Outputs/`.

Count reconciliation does not establish biological validity of the exclusions. The large excluded ribosomal clusters also contain substantial T-cell-marker expression; a completed restoration sensitivity and alternative-resource comparison are documented below. They preserve the broad non-classical-monocyte receiving pattern but change individual interaction results and do not validate blanket low-quality exclusions. Age cohort and technical batch are confounded, and inferred communication scores are relative ranks rather than direct measurements of signaling. A clean end-to-end reproduction of the current packaged code has not been independently verified.

## Annotation and resource sensitivity analyses

See [the completed sensitivity report](Results/annotation_resource_sensitivity/README.md) for all four scenarios, restored-cell counts, resource-evidence checks, and interpretation limits. Original primary outputs are preserved.

With the external archive and the CellTypist `Immune_All_Low.pkl` model available, run from the repository root:

```powershell
python Notebooks/14_ANNOTATION_RESOURCE_SENSITIVITY.py --archive-root "C:/path/to/GSE135779_transfer" --phase annotate --model "C:/path/to/Immune_All_Low.pkl"
python Notebooks/14_ANNOTATION_RESOURCE_SENSITIVITY.py --archive-root "C:/path/to/GSE135779_transfer" --phase infer --magnitude-only
python Notebooks/15_SUMMARIZE_ANNOTATION_RESOURCE_SENSITIVITY.py
```

Script 14 uses the cohort files documented for script 13. Cached predictions and donor scores are stored under ignored `_local_submission/sensitivity_cache/`; compact research outputs are saved under `Results/annotation_resource_sensitivity/`. Existing donor score files are skipped; changed inputs or settings require a fresh cache. The recorded magnitude-only equivalence check is explained in the sensitivity protocol. Script 15 requires the original saved model results and donor metadata and verifies baseline reproduction before exporting comparisons.

## Current primary annotation revision

After the sensitivity workflow above:

```powershell
python Notebooks/16_REVIEW_RESTORED_T_CELLS.py --archive-root "C:/path/to/GSE135779_transfer"
python Notebooks/17_BUILD_REVISED_PRIMARY.py --archive-root "C:/path/to/GSE135779_transfer" --publish
```

Script 16 reviews cells without loading communication outcomes. Script 17 uses the reviewed restored annotation set as primary, refits primary and female-only models, checks selected HC3 fits independently, and regenerates figures and tables. The revision is retrospective. Exact labels, input hashes, remaining uncertainty, and publication-copy behavior are described in [the current primary package](Results/revised_primary/README.md).

## License

Code and project documentation are provided under the [MIT License](LICENSE). Third-party datasets and publications remain subject to their original terms.
