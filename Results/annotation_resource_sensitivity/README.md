# Annotation and interaction-resource sensitivity results

**Version note:** this report retains the original-exclusion analysis as its comparison baseline. Following the additional [cell review](../annotation_review/README.md), the restoration/consensus scenario was adopted retrospectively as the [current primary analysis](../revised_primary/README.md). Use that folder for current-primary comparisons, female-only results, and manuscript numbers.

Both concerns were investigated using the recovered expression matrices and all 56 donors. Non-classical monocytes remain the most frequent receiving cell type among significant interactions in every scenario. Individual interaction findings depend on annotation and resource choices. These results do **not** validate blanket exclusion of the original ribosomal clusters as low-quality.

## Completed comparisons

| Scenario | Tested interactions | FDR < 0.05 | Original significant interactions still significant | Significant interactions targeting non-classical monocytes |
| --- | ---: | ---: | ---: | ---: |
| Original annotations, consensus | 2,527 | 123 | 123/123 | 54 |
| Restored T-cell subset, consensus | 2,374 | 171 | 104/119 testable | 57 |
| Original annotations, bundled CellPhoneDB | 531 | 36 | 10/14 testable | 15 |
| Restored T-cell subset, bundled CellPhoneDB | 514 | 29 | 9/14 testable | 15 |

Each scenario repeats donor-level LIANA inference and applies BH correction across its own eligible tests. Different resources change the rank universe and available pairs; counts of significant tests are therefore not directly comparable measures of signaling strength.

## Cluster exclusions

CellTypist independently classified all 348,274 post-QC cells without majority voting. A targeted sensitivity restored excluded cells only if they expressed both CD3D and CD3E, had a clearly CD4/helper or CD8/cytotoxic reference class, and had maximum model probability >=0.5. Original retained labels and RBC/platelet exclusions were preserved. See [the protocol](PROTOCOL.md) and script 14 for the explicit reference-class mapping.

This restored 76,480 pediatric and 17,434 adult cells, increasing the retained total from 204,696 to 298,610. The four excluded clusters had median ribosomal count fractions of approximately 46–50%; retained pediatric CD4 T cells had a median of approximately 43%. Fractions use genes available in each cohort matrix, not the complete original raw gene universe. Reference predictions and high ribosomal expression do not independently establish cell identity or cell quality.

Annotation sensitivity coefficients correlated with baseline at r = 0.910 across 2,371 shared tests. Of the original 123 significant results, 119 remained testable, 104 remained significant, 15 lost significance, and four were no longer testable. Two of the 119 changed coefficient sign. These checks support the broad receiving-cell pattern but show that the exact original significant set is not invariant to exclusions. They examine a defined subset of plausible T cells, not all possible revised annotations.

## Unusual resource pairs

The installed consensus resource includes both flagged pairs. Their entries were traced in the [CellTalkDB source database](https://github.com/ZJUFanLab/CellTalkDB/tree/master/database); source hashes and extracted evidence identifiers are recorded in `flagged_pair_evidence.json`.

- **ARF6–SMAP1:** the listed evidence, [PMID 15659652](https://pubmed.ncbi.nlm.nih.gov/15659652/), describes intracellular ARF6 regulation and clathrin-dependent endocytosis. It does not establish extracellular ligand–receptor communication between cells.
- **GNAI2–S1PR4:** the listed evidence, [PMID 32196115](https://pubmed.ncbi.nlm.nih.gov/32196115/), is the SingleCellSignalR methods/resource paper. That citation alone does not establish primary experimental evidence for extracellular communication by this pair. This check does not exhaust every underlying source of every database containing the pair.

Neither pair is in the CellPhoneDB resource bundled with LIANA 1.8.1. This alternative resource contains 1,223 pairs versus 4,624 in consensus. Only 14 of the original 123 significant source–target/pair combinations are represented with matching complex names and orientation; the other 109 are **untested in this comparison**, not demonstrated false. Resource-only coefficient correlation is r = 0.860 across 342 shared tests; combined sensitivity is r = 0.775 across 331. This is one established-resource sensitivity, not exhaustive manual validation of extracellular evidence for every pair or a test of the latest external CellPhoneDB release.

The flagged entries should be reported as resource-derived associations, not established intercellular signaling mechanisms. Their original rows are retained transparently, with the alternative-resource analysis reported alongside them.

## Verification and files

All 56 donors' baseline interaction sets and magnitude ranks were also checked against archived donor outputs. The maximum rank difference was 3.39e-21 (`donor_baseline_reproduction.csv`).

The complete rerun of original annotations plus consensus reproduced all 2,527 baseline models: maximum coefficient difference approximately 1.1e-15 and maximum BH q-value difference approximately 3.9e-12 (`baseline_reproduction.json`). A first-donor comparison in all four scenarios showed identical magnitude ranks using 1,000 permutations and the supported magnitude-only mode (`magnitude_only_equivalence.csv`). Remaining inference omitted unused permutation specificity scores; the outcome remains `1 - magnitude_rank`.

- `scenario_summary.csv`: the comparison above, including directional agreement and donor ranges.
- `*_models.csv`, `*_matched_comparison.csv`, `*_target_summary.csv`: full scenario results and matched comparisons.
- `excluded_cluster_review_summary.csv`, `restoration_by_group.csv`, and cohort QC/reference-label files: annotation evidence and restoration counts.
- `*_restored_cell_ids.csv.gz`: exact restored cells.
- `classifier_provenance.json`, `resource_provenance.json`, `analysis_provenance.json`: model/resource/code provenance.
- `sensitivity_coefficient_comparison.png`: baseline versus each sensitivity on shared tests; blue points met baseline FDR <0.05, gray points did not; diagonal lines mark equal coefficients.

The numerical results in this historical comparison remain unchanged. The root manuscript and visible figures/tables now use the reviewed restoration/consensus analysis as primary, with the original exclusions as sensitivity. The evidence does not establish that every remaining excluded cell was unusable. Batch confounding and lack of independent biological validation remain unresolved.
