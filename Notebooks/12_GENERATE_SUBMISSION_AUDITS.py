"""Export submission evidence from saved results, without refitting any model.

Optional --child-adata/--adult-adata accept original pre-exclusion objects. Exact
cluster counts cannot be reconstructed from retained-cell summaries alone.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]
QC = BASE / "Results" / "qc"
TABLES = BASE / "Results" / "manuscript_tables"
SEVERITY = BASE / "Results" / "severity_analysis"
CLINICAL_FIELDS = ["Race", "Ethnicity", "MMF", "OS", "MTX", "Plaquenil", "Collection_year"]


def export_clinical_covariates():
    raw = pd.read_excel(BASE / "suppdata.xlsx", sheet_name="ST1b-Clinical information ", header=None)
    header = raw.index[raw.iloc[:, 0].eq("Names")].item()
    clinical = pd.read_excel(BASE / "suppdata.xlsx", sheet_name="ST1b-Clinical information ", header=header)
    metadata = pd.read_csv(SEVERITY / "patient_clinical_metadata.csv")
    clinical = clinical[["Names"] + CLINICAL_FIELDS].rename(columns={"Names": "patient_name"})
    merged = metadata.drop(columns=CLINICAL_FIELDS, errors="ignore").merge(
        clinical, on="patient_name", how="left", validate="one_to_one", indicator=True
    )
    if not merged["_merge"].eq("both").all():
        raise ValueError("Not all analyzed donors match the source clinical sheet")
    merged = merged.drop(columns="_merge")
    merged.to_csv(SEVERITY / "patient_clinical_metadata.csv", index=False)
    audit = merged[["sample", "patient_name", "Groups"] + CLINICAL_FIELDS]
    audit.to_csv(QC / "patient_clinical_covariates.csv", index=False)
    records = []
    for group, donors in audit.groupby("Groups", sort=True):
        for field in CLINICAL_FIELDS:
            values = donors[field].astype("string").fillna("missing")
            # Preserve original codes. In particular, ND is not a negative result.
            for level, count in values.value_counts(dropna=False).items():
                missing = str(level).strip().lower() in {"nd", "na", "not reported", "missing", ""}
                records.append({"group": group, "field": field, "source_value": level,
                                "missing_or_not_reported": missing, "n_donors": count,
                                "group_total": len(donors)})
    pd.DataFrame(records).to_csv(QC / "clinical_covariate_summary.csv", index=False)


def export_model_coverage():
    primary = pd.read_csv(SEVERITY / "age_diseasestatus_interaction.csv")
    primary.groupby("n_observations").size().rename("n_interactions").reset_index().to_csv(
        QC / "primary_model_donor_coverage.csv", index=False
    )
    metadata = pd.read_csv(SEVERITY / "patient_clinical_metadata.csv")
    counts = pd.concat([pd.read_csv(QC / f"{c}_sample_celltype_counts.csv") for c in ("child", "adult")])
    if counts.duplicated(["sample", "cell_type"]).any():
        raise ValueError("Duplicate donor-cell-type counts")
    cell_types = sorted(counts["cell_type"].unique())
    grid = pd.MultiIndex.from_product([metadata["sample"], cell_types], names=["sample", "cell_type"]).to_frame(index=False)
    grid = grid.merge(counts[["sample", "cell_type", "n_cells"]], how="left", validate="one_to_one")
    grid["n_cells"] = grid["n_cells"].fillna(0).astype(int)
    grid["eligible_at_10_cells"] = grid["n_cells"] >= 10
    grid = grid.merge(metadata[["sample", "Groups"]], on="sample", validate="many_to_one")
    grid.to_csv(TABLES / "supplementary_donor_cell_type_coverage.csv", index=False)

    # The source contains counts from the original donor-score analysis.
    means = pd.read_csv(TABLES / "significant_interaction_group_means.csv")
    significant = pd.read_csv(SEVERITY / "primary_significant_interactions_fdr05.csv")
    count_cols = ["n_cH", "n_cSLE", "n_aH", "n_aSLE"]
    enriched = significant.merge(means[["interaction_id"] + count_cols], on="interaction_id", validate="one_to_one", how="left")
    if enriched[count_cols].isna().any().any() or not enriched[count_cols].sum(axis=1).eq(enriched["n_observations"]).all():
        raise ValueError("Saved group counts do not reconcile to model observations")
    if not enriched[count_cols].ge(5).all().all():
        raise ValueError("Significant result violates the donor eligibility rule")
    enriched["direction"] = enriched["interaction_coef"].map(
        lambda x: "positive (more positive pediatric SLE-associated change)" if x > 0
        else "negative (more negative pediatric SLE-associated change)"
    )
    enriched.to_csv(TABLES / "supplementary_all_fdr05_interactions.csv", index=False)


def export_exclusion_inventory(paths):
    records = []
    totals = pd.read_csv(QC / "annotation_exclusion_summary.csv").set_index("cohort")
    for cohort, path in paths.items():
        mapping = pd.read_csv(QC / f"{cohort}_cluster_to_celltype_map.csv", dtype={"cluster": str})
        bad = {"Ribosomal/Low-quality", "Mitochondrial-high/Low-quality", "Mixed platelet-monocyte/Low-quality", "RBCs", "Platelets"}
        evidence_path = QC / f"{cohort}_annotation_cluster_counts.csv"
        if path.exists():
            import scanpy as sc
            from qc_utils import export_annotation_evidence
            adata = sc.read_h5ad(path, backed="r")
            try:
                adata.obs["leiden"] = adata.obs["leiden_res_0.8"].astype(str)
                export_annotation_evidence(adata, dict(zip(mapping["cluster"], mapping["cell_type"])), bad, cohort, BASE)
            finally:
                adata.file.close()
        counts = pd.read_csv(evidence_path, dtype={"cluster": str}) if evidence_path.exists() else None
        if counts is not None:
            excluded_total = counts.loc[counts["excluded"], "n_cells"].sum()
            if excluded_total != totals.loc[cohort, "cells_excluded"]:
                raise ValueError(f"{cohort}: excluded clusters do not reconcile to saved total")
        for row in mapping.itertuples(index=False):
            if row.cell_type not in bad:
                continue
            matches = counts.loc[counts["cluster"].eq(row.cluster), "n_cells"] if counts is not None else []
            records.append({"cohort": cohort, "cluster": row.cluster, "assigned_cell_type": row.cell_type,
                            "n_cells": int(matches.iloc[0]) if len(matches) == 1 else pd.NA,
                            "count_status": "verified from pre-exclusion object" if len(matches) == 1 else "unavailable: pre-exclusion object required",
                            "cohort_exclusions_total": int(totals.loc[cohort, "cells_excluded"]),
                            "marker_table_available": (QC / f"{cohort}_marker_genes_all_clusters.csv.gz").exists()})
    pd.DataFrame(records).to_csv(QC / "annotation_exclusion_evidence_inventory.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child-adata", type=Path, default=BASE / "child_individual_h5ad" / "adata_child_processed_leiden08.h5ad")
    parser.add_argument("--adult-adata", type=Path, default=BASE / "adult_individual_h5ad" / "adata_adult_processed_selected_leiden.h5ad")
    args = parser.parse_args()
    export_clinical_covariates()
    export_model_coverage()
    export_exclusion_inventory({"child": args.child_adata, "adult": args.adult_adata})
    print("Updated clinical, model-coverage, and exclusion-evidence audits; no models refitted.")


if __name__ == "__main__":
    main()
