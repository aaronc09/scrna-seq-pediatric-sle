"""Read archived AnnData files without loading their full expression matrices.

Reconcile observations to the manuscript QC outputs and summarize canonical
markers in bounded CSR chunks. This audit does not relabel cells or refit models.
Requires pandas, numpy, and h5py. All source files are opened read-only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
QC = BASE / "Results" / "qc"
MARKERS = ["CD3D", "CD3E", "TRAC", "IL7R", "CCR7", "TCF7", "CD4", "CD8A", "CD8B",
           "MS4A1", "CD79A", "CD37", "NKG7", "GNLY", "KLRD1", "LYZ", "LST1", "FCGR3A",
           "S100A8", "S100A9", "FCER1A", "CD1C", "GZMB", "TCF4", "PLD4", "IGJ", "MZB1",
           "IGLL5", "ISG15", "IFIT1", "HBA1", "HBA2", "HBB", "ALAS2", "PPBP", "PF4",
           "TUBB1", "LTB", "JUNB"]
BAD = {"Ribosomal/Low-quality", "Mitochondrial-high/Low-quality",
       "Mixed platelet-monocyte/Low-quality", "RBCs", "Platelets"}


def decode(node):
    if isinstance(node, h5py.Group):
        return pd.Categorical.from_codes(node["codes"][:], decode(node["categories"]))
    return node.asstr()[:] if node.dtype.kind in ("O", "S") else node[:]


def check_cohort(root, cohort):
    name = "adata_child_processed_leiden08.h5ad" if cohort == "child" else "adata_adult_processed_selected_leiden.h5ad"
    notebook = "01B_GSE135779_CHILD_CELLTYPE_LABELLING" if cohort == "child" else "01E_GSE135779_ADULT_CELLTYPE_LABELLING"
    source = root / f"{cohort}_individual_h5ad" / name
    marker_path = root / "Notebook_Outputs" / notebook / "marker_genes_all_clusters.csv"
    mapping = pd.read_csv(QC / f"{cohort}_cluster_to_celltype_map.csv", dtype={"cluster": str}).set_index("cluster")["cell_type"]
    with h5py.File(source, "r") as handle:
        group = handle["obs"]
        cols = ["sample", "leiden_res_0.8", "n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_score"]
        obs = pd.DataFrame({k: decode(group[k]) for k in cols}, index=decode(group[group.attrs["_index"]]))
        obs.index.name = "cell_id"
        obs["sample"] = obs["sample"].astype(str)
        obs["cluster"] = obs["leiden_res_0.8"].astype(str)
        obs["assigned_cell_type"] = obs["cluster"].map(mapping)
        obs["excluded"] = obs["assigned_cell_type"].isin(BAD)
        if obs["assigned_cell_type"].isna().any() or not obs.index.is_unique:
            raise ValueError("Unmapped clusters or nonunique cell IDs")
        post = pd.read_csv(QC / f"{cohort}_postfilter_sample_qc.csv").set_index("sample")["n_cells"]
        observed = obs.groupby("sample", observed=True).size()
        pd.testing.assert_series_equal(post.sort_index(), observed.sort_index(), check_names=False, check_index_type=False)
        counts = obs.loc[~obs["excluded"]].groupby(["sample", "assigned_cell_type"], observed=True).size()
        saved = pd.read_csv(QC / f"{cohort}_sample_celltype_counts.csv").set_index(["sample", "cell_type"])["n_cells"]
        counts.index = counts.index.set_names(saved.index.names)
        pd.testing.assert_series_equal(counts.sort_index(), saved.sort_index(), check_names=False, check_index_type=False)
        audit = pd.read_csv(QC / "annotation_exclusion_by_sample.csv")
        expected = audit.loc[audit["cohort"].eq(cohort)].set_index("sample")["excluded_by_annotation"]
        actual = obs.groupby("sample", observed=True)["excluded"].sum()
        pd.testing.assert_series_equal(expected.sort_index(), actual.sort_index(), check_names=False, check_index_type=False)

        keys = ["cluster", "assigned_cell_type", "excluded"]
        clusters = obs.groupby(keys, observed=True).size().rename("n_cells").reset_index()
        for key in cols[2:]:
            clusters[f"median_{key}"] = clusters["cluster"].map(obs.groupby("cluster", observed=True)[key].median())
        clusters = clusters.sort_values("cluster", key=lambda x: x.astype(int))
        clusters.to_csv(QC / f"{cohort}_annotation_cluster_counts.csv", index=False)
        obs.groupby(["sample"] + keys, observed=True).size().rename("n_cells").reset_index().to_csv(
            QC / f"{cohort}_annotation_counts_by_sample_cluster.csv", index=False)
        obs.loc[obs["excluded"], ["sample"] + keys].reset_index().to_csv(QC / f"{cohort}_excluded_cell_ids.csv.gz", index=False)

        genes = decode(handle["var"][handle["var"].attrs["_index"]])
        marker_table = pd.read_csv(marker_path, dtype={"group": str}, float_precision="round_trip")
        if set(marker_table["group"]) != set(obs["cluster"]):
            raise ValueError("Marker table clusters do not match archived object")
        for _, rows in marker_table.groupby("group"):
            if set(rows["names"]) != set(genes) or len(rows) != len(genes):
                raise ValueError("Full marker ranking gene universe does not match archived object")
        marker_table.to_csv(QC / f"{cohort}_marker_genes_all_clusters.csv.gz", index=False)
        marker_table.groupby("group", sort=False).head(20).to_csv(QC / f"{cohort}_top20_cluster_markers.csv", index=False)

        selected = [g for g in MARKERS if g in set(genes)]
        lookup = np.full(len(genes), -1, dtype=int)
        positions = {g: i for i, g in enumerate(genes)}
        for i, gene in enumerate(selected):
            lookup[positions[gene]] = i
        labels = sorted(obs["cluster"].unique(), key=int)
        codes = pd.Categorical(obs["cluster"], categories=labels).codes
        n_cells = np.bincount(codes)
        sums = np.zeros((len(labels), len(selected)), dtype=np.float64)
        detected = np.zeros_like(sums, dtype=np.int64)
        t_both = np.zeros(len(labels), dtype=np.int64)
        matrix = handle["X"]
        if matrix.attrs["encoding-type"] != "csr_matrix":
            raise ValueError("Expected CSR expression matrix")
        ptr = matrix["indptr"][:]
        for start in range(0, len(obs), 4000):
            stop = min(start + 4000, len(obs))
            left, right = int(ptr[start]), int(ptr[stop])
            index = matrix["indices"][left:right]
            local_gene = lookup[index]
            keep = local_gene >= 0
            entries = np.flatnonzero(keep)
            local_row = np.searchsorted(ptr[start:stop+1], left + entries, side="right") - 1
            values = matrix["data"][left:right][keep]
            block = np.zeros((stop-start, len(selected)), dtype=np.float32)
            block[local_row, local_gene[keep]] = values
            for code in np.unique(codes[start:stop]):
                subset = block[codes[start:stop] == code]
                sums[code] += subset.sum(axis=0, dtype=np.float64)
                detected[code] += (subset > 0).sum(axis=0)
                t_both[code] += ((subset[:, selected.index("CD3D")] > 0) & (subset[:, selected.index("CD3E")] > 0)).sum()
        records = [{"cluster": label, "assigned_cell_type": mapping[label], "marker": gene,
                    "n_cells": int(n_cells[i]), "mean_expression": sums[i,j]/n_cells[i],
                    "fraction_expressing": detected[i,j]/n_cells[i]}
                   for i,label in enumerate(labels) for j,gene in enumerate(selected)]
        expression = pd.DataFrame(records)
        expression.to_csv(QC / f"{cohort}_all_cluster_marker_expression.csv", index=False)
        clusters["fraction_CD3D_and_CD3E"] = [t_both[labels.index(c)] / n_cells[labels.index(c)] for c in clusters["cluster"]]
        clusters.to_csv(QC / f"{cohort}_annotation_cluster_counts.csv", index=False)
        # Independent expression-level cross-check against the manuscript's saved
        # retained-lineage summaries. Multiple clusters are weighted by cell count.
        old = pd.read_csv(QC / f"{cohort}_canonical_marker_validation.csv")
        deltas = []
        for row in old.loc[old["present"]].itertuples(index=False):
            subset = expression.loc[expression.assigned_cell_type.eq(row.cell_type) & expression.marker.eq(row.marker)]
            if subset.empty:
                continue
            mean = np.average(subset.mean_expression, weights=subset.n_cells)
            frac = np.average(subset.fraction_expressing, weights=subset.n_cells)
            deltas.append({"cell_type": row.cell_type, "marker": row.marker,
                           "absolute_mean_difference": abs(mean-row.mean_expression),
                           "absolute_fraction_difference": abs(frac-row.fraction_expressing)})
        pd.DataFrame(deltas).to_csv(QC / f"{cohort}_archive_expression_crosscheck.csv", index=False)
        max_mean = max(x["absolute_mean_difference"] for x in deltas)
        max_fraction = max(x["absolute_fraction_difference"] for x in deltas)
        if max_mean > 1e-5 or max_fraction > 1e-10:
            raise ValueError("Archived expression does not reproduce saved canonical-marker summaries")
        result = {"source_h5ad": source.relative_to(root).as_posix(), "source_bytes": source.stat().st_size,
                  "source_mtime_utc": datetime.fromtimestamp(source.stat().st_mtime, timezone.utc).isoformat(),
                  "marker_source": marker_path.relative_to(root).as_posix(), "marker_sha256": hashlib.sha256(marker_path.read_bytes()).hexdigest(),
                  "cells_before_annotation": len(obs), "cells_excluded": int(obs.excluded.sum()),
                  "cells_retained": int((~obs.excluded).sum()), "donors": len(post),
                  "all_donor_counts_match": True, "all_retained_donor_celltype_counts_match": True,
                  "all_donor_exclusion_counts_match": True, "marker_gene_universe_matches": True,
                  "max_retained_marker_mean_difference": max_mean,
                  "max_retained_marker_detection_difference": max_fraction,
                  "biological_exclusion_validity": "not established by count reconciliation"}
        print(cohort, json.dumps(result), flush=True)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-root", type=Path, required=True)
    args = parser.parse_args()
    results = {c: check_cohort(args.archive_root, c) for c in ("child", "adult")}
    results["source_path_base"] = "Paths are relative to the external archive root; source matrices are not bundled."
    results["audit_utc"] = datetime.now(timezone.utc).isoformat()
    (BASE / "Results" / "provenance" / "archived_annotation_audit.json").write_text(json.dumps(results, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
