CORRECTED REVIEWER REVISION - 2026-10-06

Authoritative current inference and model outputs are in this folder.
The uncorrected submitted analyses in Results/revised_primary/ and
Results/annotation_resource_sensitivity/ are retained for comparison.
Their cell labels, cell counts, annotation evidence, UMAP and cohort table
remain relevant; their LIANA magnitude-based model results are historical.

Corrected primary: 220 significant / 2335 estimable tests.
Old significant results retained: 116 / 171.
All original 2,374 tests retain donor-availability eligibility; 39 corrected
score outcomes are constant and have no estimable test.

READ FIRST
analysis_protocol.txt: revision-stage protocol and implementation clarifications.
summary.json: machine-readable numerical summary.
primary_models.csv: authoritative fresh-inference primary estimates.
supplementary_all_significant.csv: corrected FDR < 0.05 combinations (Supplemental Table S3A).
annotation_shared_significant.csv: significant under both annotation approaches.
software_versions.json and provenance.json: exact versions and input/code hashes.
inference_verification.csv: all 56 fresh primary and return-all inference checks.
independent_model_verification.csv: independent statsmodels HC3 checks.

FIGURES AND TABLES
figure_table_numbering.txt: main and supplemental numbering and panel-to-file map.
Figures 1-5 are in Results/manuscript_figures/ (PNG and available SVG formats).
figure2_primary_analysis.svg/png and figure3_female_only_comparison.svg/png are Figures 2 and 3.
table2_editable.html: native selectable Table 2; copy into the journal document.
table2_panel_A.csv and table2_panel_B.csv: the same native table values.
Results/manuscript_tables/table2_interaction_estimates.png: both corrected panels, 18-point body text, 600 dpi.
Regenerate this PNG alone with: python Notebooks/render_table2_png.py
table2_corrected.csv: full corrected/submitted example comparison.
Results/manuscript_figures/figure4_individual_donor_scores.svg/png: Figure 4.
Results/manuscript_figures/figure5_pediatric_SLE_by_batch.svg/png: Figure 5.
Original Figure 1 and Table 1 remain unchanged; labels and coordinates were fixed.
PNG exports are 400 dpi; SVG provides scalable vector output.

DONOR AVAILABILITY
all_primary_tests_donor_availability.csv.gz: 2,374 original tests x 56 donors.
significant_donor_missingness_reasons.csv.gz: old/new significant union.
donor_availability_by_group.csv: per-group counts for the entire original family.
Flags retain overlapping population/expression failures. Population inadequacy
takes precedence in the human-readable missing_reason field.
two_part_models.csv: expression detectability and conditional-score coefficients
and separately adjusted q values. Inadequate populations are not scored as zero.
return_all_original_family_models.csv: return-all sensitivity on original family.
common_support_fixed_universe.csv: group-blind fixed support universe.
common_support_models.csv: estimable common-support tests.

OTHER SENSITIVITIES
batch_adjusted_models.csv: estimability diagnostics and additive batch models.
group_by_batch_counts.csv and table2_pediatric_SLE_batch_tests.csv.
table2_individual_donor_scores.csv and table2_leave_one_donor_out.csv.
female_models.csv, baseline_models.csv (original labels/consensus),
resource_models.csv (original labels/CellPhoneDB), and combined_models.csv
(revised labels/CellPhoneDB); *_comparison.csv files compare with primary.
annotation_models.csv is the cached-score reconstruction check; use
primary_models.csv from fresh inference for reported primary estimates.
old_vs_corrected_comparison.csv and previously_tested_now_nonestimable.csv.
submitted_104_after_correction.csv tracks the historical annotation overlap.

REPRODUCTION
Use the recorded Python environment, including LIANA 1.8.1. Script 18 patches
ranking only in its own process; running unpatched earlier scripts by themselves
will reproduce the historical bug. Do not broadly upgrade packages to isolate
this correction.

Prerequisites: earlier project stages through scripts 14-17, the two archived
post-QC expression objects, local sensitivity_cache/*_cell_review.csv.gz, and
all four scenarios of saved per-method donor scores. The repository does not
include the large local expression/cache files. Regenerate them with the earlier
pipeline or supply the original local archive; the compact outputs here support
inspection but do not replace raw input acquisition.

python Notebooks/18_REVIEWER_REVISION.py cached
python Notebooks/18_REVIEWER_REVISION.py infer --archive-root "C:/path/to/GSE135779_transfer"
python Notebooks/19_SUMMARIZE_REVIEWER_REVISION.py

Script 18 verifies old cached scores, fixes each unique magnitude-score rank
once, freezes checked resource snapshots, and reruns donor inference from X.
The inference stage checkpoints each donor only after primary, availability,
and return-all checks complete. Fresh primary output is authoritative because
CSV round trips can split tied raw scores at very small precision differences.
The other annotation/resource scenarios use corrected aggregation of unchanged
per-method scores. Script 19 refits models, checks them, and creates this package's
scientific outputs. Existing manuscript and response prose are revision artifacts;
the scripts do not automatically rewrite them on future runs.

INTERPRETATION AND HANDOFF
The scores are relative rank quantities, not absolute signaling intensity.
The revision does not eliminate age-batch confounding, small-sample uncertainty,
annotation selection, uncertain resource pairs, or the need for replication.
reference_verification.txt records the 29-reference review.
replication_search.txt records the limited independent-cohort search.
The local manuscript_draft.txt and response_to_reviewers.txt contain the revision.
This task did not upload changes to GitHub; publish the revised code/results
before making a public availability claim in the resubmitted manuscript.
