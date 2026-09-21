# Annotation and interaction-resource sensitivity protocol

This is a retrospective robustness analysis, not a preregistered analysis. The
rules below were selected before inspecting sensitivity-model outcomes. Original
annotations, expression files, primary results, figures, and tables are preserved.

## Annotation comparison

Apply the locally available CellTypist `Immune_All_Low.pkl` model to all post-QC
cells in both cohorts using independent per-cell predictions (no majority voting).
Review predicted labels and cell-level QC, including ribosomal count fractions.
These fractions use the available cohort gene intersection as the denominator,
not every gene in the original raw library. Predictions are supporting evidence,
not experimentally established cell identities or calibrated quality scores.

For a targeted restoration sensitivity, retain every originally retained cell
with its original label. From the four excluded ribosomal clusters (child 0/1,
adult 0/4), restore cells only when both CD3D and CD3E are detected and a clearly
helper/CD4 or cytotoxic/CD8 reference class has a maximum model probability of at
least 0.5. Map these classes to the original broad CD4/CD8 labels, using the
explicit class lists in notebook script 14. Do not force ambiguous T, NKT, MAIT,
gamma-delta, or thymocyte predictions into these labels. Keep RBC/platelet
exclusions unchanged. All cells already passed the original preprocessing QC;
high ribosomal expression alone is not an additional exclusion rule. Report
probability >=0.8 counts as a diagnostic without choosing a cutoff based on the
communication results. This sensitivity tests a defined subset of plausible
excluded T cells, not every possible reannotation of every lineage.

## Resource comparison

Use the CellPhoneDB resource bundled with LIANA 1.8.1 as an alternative resource
for every donor, irrespective of statistical significance. Keep the inference
method (LIANA rank aggregation), expression/cell thresholds, permutations, and
model unchanged. This is an established-resource comparison, not a claim that
every entry has been manually revalidated here. Freeze both resource tables and
the installed bundle hash. It is not the latest external CellPhoneDB release.

## Four comparisons

1. Original annotations + consensus resource (reproduction control).
2. Restored T-cell subset + consensus resource (annotation sensitivity).
3. Original annotations + bundled CellPhoneDB (resource sensitivity).
4. Restored T-cell subset + bundled CellPhoneDB (combined sensitivity).

Run donor-level inference afresh, including rank aggregation within each resource.
Fit the original disease-by-age-group OLS model with HC3 errors, requiring at
least five donors per group, and recompute BH FDR across each scenario's eligible
tests. Compare coefficients on matched interactions, directions, significant-set
overlap, donor counts, and non-classical-monocyte targets. The scenarios use
different rank universes, so score differences are not absolute signaling changes.
No selective removal of significant rows or reuse of their old q values is used.

Implementation note: the first donor was run with the original 1,000 permutations
in all four scenarios. LIANA's magnitude-only mode (`n_perms=None`) then reproduced
every magnitude rank exactly in each scenario (see `magnitude_only_equivalence.csv`).
The remaining inference skips unused permutation-based specificity scores using
that supported mode; the downstream outcome is exclusively `1 - magnitude_rank`.
Every donor's original-annotation consensus scores are rerun, and the complete
baseline models must match the saved original analysis before reporting comparisons.

Intermediate cell-level predictions and donor scores are kept under the ignored
`_local_submission/sensitivity_cache/`; compact summaries and restored-cell IDs
are saved here. Source AnnData objects are opened read-only. Resource checks and
CellTypist cannot resolve all annotation uncertainty or batch confounding.
