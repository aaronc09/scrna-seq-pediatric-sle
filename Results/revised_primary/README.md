# Current primary analysis

This folder is the authoritative numerical source for the revised root manuscript. It uses the reviewed conservative T-cell restoration and LIANA 1.8.1 consensus resource. Original full-cluster exclusions are now a sensitivity comparison. This is a retrospective annotation revision, not a preregistered or externally validated annotation.

| Result | Revised primary |
| --- | ---: |
| Donors | 56 |
| Pediatric cells | 233,770 |
| Adult cells | 64,840 |
| Total retained cells | 298,610 |
| Eligible interaction tests | 2,374 |
| FDR <0.05 | 171 |
| Positive / negative coefficients at FDR <0.05 | 37 / 134 |
| FDR <0.10 | 239 |
| Significant interactions targeting non-classical monocytes | 57 |
| Female-only tests / FDR <0.05 | 2,373 / 145 |
| Primary–female coefficient correlation | 0.994839 |
| Primary significant results retained in female-only analysis | 141 / 171 |

The [annotation review](../annotation_review/README.md) explains the supporting evidence and unresolved uncertainty. The selection restores 76,480 pediatric and 17,434 adult cells; 45,420 cells from the original T-rich clusters remain unresolved and excluded. Their exclusion is not justified by claiming they are all low-quality. The reference-confidence and marker-detection requirements can still introduce selection bias.

## Files and verification

- `primary_models.csv`: all revised primary tests, confidence intervals, p/q values, group means, and contributing-donor counts.
- `female_only_models.csv`, `female_only_comparison.csv`: newly fitted female-only models and matched comparison.
- `*_primary_cell_labels.csv.gz`: exact cell identifiers, original and revised labels, restoration flag, and inclusion status.
- `qc/annotation_accounting.csv`: complete donor-level post-QC, restoration, exclusion, and retention accounting.
- `sensitivity_summary.csv` and `*_comparison.csv`: comparisons using revised-primary denominators. These differ from the historical original-relative comparisons under `annotation_resource_sensitivity/`.
- `../manuscript_tables/`: the single current set of main tables and supplementary data. The original seven Table 2 example combinations were retained and all remain significant; examples were not replaced to maximize significance.
- `../manuscript_figures/`: the single current research figure set, Figures 1–3 as PNG and JPEG, each exported at 600 DPI. Figure 1 uses archived pre-exclusion UMAP coordinates and revised labels, without fitting a new embedding.
- `summary.json`, `provenance.json`, `input_sha256.json`: machine-readable results, methods, and frozen input hashes.
- `independent_HC3_verification.csv`: selected primary and female-only fits checked against statsmodels formula-based HC3 estimates. Full primary coefficients and q values also reproduce the previously completed annotation scenario.

The exact annotation rule matches the completed restoration scenario, so all 56 donor inferences were reused rather than changed. Primary models were refit and checked; the female-only models, tables, figures, accounting, and primary-relative comparisons were generated for this revision. The prior baseline rerun had already reproduced archived donor inference and original statistical results.

With revised annotations, the bundled CellPhoneDB sensitivity identifies 29 significant interactions among 514 tests. Only 17 of the revised primary 171 significant interactions are testable in that comparison, and 10 remain significant; 154 are untested, not disproven. Non-classical monocytes remain the most frequent significant receiving population. The unusual GNAI2–S1PR4 and ARF6–SMAP1 entries remain qualified as resource-derived associations, not established intercellular mechanisms.

## Reproduction

After scripts 13–15 and the required local inference cache are available:

```powershell
python Notebooks/16_REVIEW_RESTORED_T_CELLS.py --archive-root "C:/path/to/GSE135779_transfer"
python Notebooks/17_BUILD_REVISED_PRIMARY.py --archive-root "C:/path/to/GSE135779_transfer" --publish
```

Script 17 writes this versioned package. `--publish` refreshes visible copies in `Results/manuscript_figures/`, `Results/manuscript_tables/`, and `Journal_Submission/Figures/`, preserving prior presentation files under ignored `_local_submission/pre_primary_revision/`. It does not edit manuscript prose. Manuscript drafts and submission packages are maintained separately from the research repository. The optional journal exports and local backups are ignored by Git.

After publishing, script 17 removes the verified duplicate figure/table staging folders from this directory and combines numbered historical run manifests in `Results/provenance/run_manifests.json`. Model outputs and original biological results are unchanged. A build without `--publish` retains staging folders until publication. The recorded generation-script hash describes the code used for the original numerical run; later housekeeping changes do not constitute a new scientific run.

Original model/QC outputs under `Results/severity_analysis/`, `Results/sensitivity_analysis/`, and `Results/qc/` are preserved as historical inputs and comparison evidence. Do not use their original-primary counts for the revised manuscript. Older plotting scripts can regenerate old presentations; run script 17 with `--publish` last when rebuilding the current package.
