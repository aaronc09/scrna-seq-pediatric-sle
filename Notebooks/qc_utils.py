"""Publication QC table exporters for AnnData objects."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import silhouette_score


CANONICAL_MARKERS = {
    "B Cells": ("CD79A", "MS4A1", "CD37"),
    "CD4 T Cells": ("CD3D", "IL7R", "LTB"),
    "CD8 T Cells": ("CD3D", "CD8A", "CD8B"),
    "NK Cells": ("NKG7", "GNLY", "KLRD1"),
    "Cytotoxic T/NK Cells": ("NKG7", "GNLY", "GZMB"),
    "Classical Monocytes": ("LYZ", "S100A8", "S100A9"),
    "Non-classical Monocytes": ("LYZ", "FCGR3A", "LST1"),
    "IFN-stimulated Monocytes": ("ISG15", "IFIT1", "IFIT3"),
    "IFN-stimulated T Cells": ("ISG15", "IFIT1", "CD3D"),
    "Plasma Cells": ("IGJ", "MZB1", "IGLL5"),
    "cDCs": ("FCER1A", "CST3", "CD1C"),
    "pDCs": ("GZMB", "TCF4", "PLD4"),
}


def export_annotation_evidence(
    adata, cell_type_map, bad_types, cohort, base_dir, marker_table=None
):
    """Audit the pre-exclusion object; never infer cluster counts by subtraction.

    Call before discarding any clusters. Marker rankings must have been computed
    on this same object and Leiden solution. No new annotation is assigned here.
    """
    output = Path(base_dir) / "Results" / "qc"
    output.mkdir(parents=True, exist_ok=True)
    obs = adata.obs.copy()
    obs["cluster"] = obs["leiden"].astype(str)
    obs["assigned_cell_type"] = obs["cluster"].map(cell_type_map)
    if obs["assigned_cell_type"].isna().any():
        raise ValueError("Unmapped Leiden clusters in annotation audit")
    if not obs.index.is_unique:
        raise ValueError("Cell identifiers must be unique for the exclusion audit")
    obs["excluded"] = obs["assigned_cell_type"].isin(bad_types)
    expected_path = output / f"{cohort}_postfilter_sample_qc.csv"
    if expected_path.exists():
        expected = pd.read_csv(expected_path).set_index("sample")["n_cells"]
        observed = obs.groupby("sample", observed=True).size()
        if set(expected.index) != set(observed.index) or not (
            expected.sort_index().to_numpy() == observed.sort_index().to_numpy()
        ).all():
            raise ValueError("Pre-exclusion object does not match saved post-QC donor counts")
    keys = ["cluster", "assigned_cell_type", "excluded"]
    clusters = obs.groupby(keys, observed=True).size().rename("n_cells").reset_index()
    clusters = clusters.sort_values("cluster", key=lambda s: s.astype(int))
    for field in ("n_genes_by_counts", "total_counts", "pct_counts_mt"):
        if field in obs:
            medians = obs.groupby("cluster", observed=True)[field].median()
            clusters[f"median_{field}"] = clusters["cluster"].map(medians)
    by_sample = obs.groupby(["sample"] + keys, observed=True).size().rename("n_cells").reset_index()
    excluded = obs.loc[obs["excluded"], ["sample"] + keys].copy()
    excluded.insert(0, "cell_id", excluded.index.astype(str))
    clusters.to_csv(output / f"{cohort}_annotation_cluster_counts.csv", index=False)
    by_sample.to_csv(output / f"{cohort}_annotation_counts_by_sample_cluster.csv", index=False)
    excluded.to_csv(output / f"{cohort}_excluded_cell_ids.csv.gz", index=False)
    if marker_table is not None:
        marker_table.to_csv(output / f"{cohort}_marker_genes_all_clusters.csv.gz", index=False)
        marker_table.groupby("group", observed=True, sort=False).head(20).to_csv(
            output / f"{cohort}_top20_cluster_markers.csv", index=False
        )
    return clusters


def export_sample_qc(adata, output_path: str | Path) -> pd.DataFrame:
    """Export post-filter QC summaries with one row per donor/sample."""
    required = {"sample", "n_genes_by_counts", "total_counts", "pct_counts_mt"}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise ValueError(f"AnnData lacks required QC columns: {sorted(missing)}")
    grouped = adata.obs.groupby("sample", observed=True)
    table = grouped.agg(
        n_cells=("sample", "size"),
        median_genes_per_cell=("n_genes_by_counts", "median"),
        median_umis_per_cell=("total_counts", "median"),
        median_pct_mitochondrial=("pct_counts_mt", "median"),
    ).reset_index()
    if "predicted_doublet" in adata.obs:
        # Retained objects should contain no predicted doublets; this documents it.
        doublets = grouped["predicted_doublet"].sum().rename("retained_predicted_doublets")
        table = table.merge(doublets.reset_index(), on="sample", how="left")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table


def export_celltype_counts(adata, output_path: str | Path) -> pd.DataFrame:
    """Export donor-by-cell-type counts and within-donor proportions."""
    required = {"sample", "cell_type"}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise ValueError(f"AnnData lacks required annotation columns: {sorted(missing)}")
    table = (
        adata.obs.groupby(["sample", "cell_type"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    table["proportion_within_sample"] = table["n_cells"] / table.groupby(
        "sample", observed=True
    )["n_cells"].transform("sum")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table


def export_marker_validation(adata, output_path: str | Path) -> pd.DataFrame:
    """Summarize canonical-marker expression for each assigned cell type."""
    rows = []
    for cell_type in sorted(adata.obs["cell_type"].dropna().unique()):
        selected = adata[adata.obs["cell_type"] == cell_type]
        for marker in CANONICAL_MARKERS.get(str(cell_type), ()):
            if marker not in adata.var_names:
                rows.append({"cell_type": cell_type, "marker": marker,
                             "present": False, "mean_expression": np.nan,
                             "fraction_expressing": np.nan})
                continue
            values = selected[:, marker].X
            values = values.toarray().ravel() if sparse.issparse(values) else np.asarray(values).ravel()
            rows.append({"cell_type": cell_type, "marker": marker,
                         "present": True, "mean_expression": float(values.mean()),
                         "fraction_expressing": float((values > 0).mean())})
    table = pd.DataFrame(rows)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table


def export_integration_diagnostics(adata, output_path: str | Path, sample_size: int = 5000) -> pd.DataFrame:
    """Compare sample mixing and cell-type conservation before/after Harmony."""
    missing = {"X_pca", "X_pca_harmony"}.difference(adata.obsm)
    if missing:
        raise ValueError(f"AnnData lacks embeddings: {sorted(missing)}")
    n = min(sample_size, adata.n_obs)
    rng = np.random.default_rng(0)
    index = rng.choice(adata.n_obs, size=n, replace=False)
    rows = []
    for name, key in [("before_harmony", "X_pca"), ("after_harmony", "X_pca_harmony")]:
        embedding = np.asarray(adata.obsm[key])[index, :30]
        for label, interpretation in [
            ("sample", "lower absolute silhouette indicates better sample mixing"),
            ("cell_type", "higher silhouette indicates better cell-type conservation"),
        ]:
            labels = adata.obs[label].astype(str).to_numpy()[index]
            value = silhouette_score(embedding, labels) if len(np.unique(labels)) > 1 else np.nan
            rows.append({"embedding": name, "label": label,
                         "silhouette_score": value, "interpretation": interpretation,
                         "n_cells_sampled": n})
    table = pd.DataFrame(rows)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table


def export_qc_threshold_sensitivity(h5ad_paths, output_path: str | Path) -> pd.DataFrame:
    """Report cell retention over a prespecified QC-threshold grid."""
    import scanpy as sc

    rows = []
    for raw_path in sorted(map(Path, h5ad_paths)):
        adata = sc.read_h5ad(raw_path)
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        for min_genes in (100, 200, 300, 500):
            for max_mito in (10.0, 15.0, 20.0, 25.0):
                retained = (
                    (adata.obs["n_genes_by_counts"] >= min_genes)
                    & (adata.obs["pct_counts_mt"] < max_mito)
                )
                rows.append({"sample": str(adata.obs["sample"].iloc[0]),
                             "min_genes": min_genes,
                             "max_pct_mitochondrial": max_mito,
                             "n_cells_total": adata.n_obs,
                             "n_cells_retained": int(retained.sum()),
                             "fraction_retained": float(retained.mean())})
    table = pd.DataFrame(rows)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table
