"""
Regenerates the Harmony batch-correction QC figures (Notebook_Outputs/batch_check_*.png).

These compare a UMAP colored by sample before vs. after Harmony batch correction, for
both cohorts. "Before" recomputes a plain (non-batch-corrected) PCA -> neighbors -> UMAP
from each cohort's post-HVG checkpoint; "after" reuses the already-computed Harmony-
corrected UMAP from the final processed file (no need to recompute Harmony itself).

Run from the Notebooks/ directory (or anywhere -- BASE_DIR is resolved relative to this
file's parent), after 1A/1D have produced their respective checkpoint files.
"""
import os
import matplotlib

matplotlib.use("Agg")

import scanpy as sc
import matplotlib.pyplot as plt

from publication_utils import (
    configure_publication_notebook,
    remove_owned_outputs,
    save_publication_figure,
)

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIGURE_DIR = configure_publication_notebook(BASE_DIR, "batch_correction_qc")


def make_before_after(before_path, after_path, out_prefix, label):
    adata_before = sc.read_h5ad(before_path)
    sc.tl.pca(adata_before, svd_solver="arpack", use_highly_variable=True, random_state=0)
    sc.pp.neighbors(adata_before, n_neighbors=15, n_pcs=30, random_state=0)
    sc.tl.umap(adata_before, random_state=0)

    sc.pl.umap(adata_before, color="sample", legend_loc=None,
               title=f"{label}: UMAP by sample (before Harmony)", show=False)
    plt.tight_layout()
    save_publication_figure(FIGURE_DIR / f"{out_prefix}_umap_by_sample.png")
    plt.close()
    del adata_before

    adata_after = sc.read_h5ad(after_path)
    sc.pl.umap(adata_after, color="sample", legend_loc=None,
               title=f"{label}: UMAP by sample (after Harmony)", show=False)
    plt.tight_layout()
    save_publication_figure(
        FIGURE_DIR / f"{out_prefix}_umap_by_sample_after_correction.png"
    )
    plt.close()
    del adata_after
    print(f"Saved {label} before/after batch-check figures")


if __name__ == "__main__":
    remove_owned_outputs(
        [
            *FIGURE_DIR.glob("*.png"),
            f"{BASE_DIR}/Notebook_Outputs/batch_check_child_umap_by_sample.png",
            f"{BASE_DIR}/Notebook_Outputs/batch_check_child_umap_by_sample_AFTER_correction.png",
            f"{BASE_DIR}/Notebook_Outputs/batch_check_adult_umap_by_sample.png",
            f"{BASE_DIR}/Notebook_Outputs/batch_check_adult_umap_by_sample_AFTER_correction.png",
        ]
    )
    make_before_after(
        before_path=f"{BASE_DIR}/adata_child_combined_log1p_hvg.h5ad",
        after_path=f"{BASE_DIR}/child_individual_h5ad/adata_child_processed_leiden08.h5ad",
        out_prefix="batch_check_child",
        label="Child",
    )
    make_before_after(
        before_path=f"{BASE_DIR}/adult_individual_h5ad/adata_adult_combined_filtered2.h5ad",
        after_path=f"{BASE_DIR}/adult_individual_h5ad/adata_adult_processed_leiden08.h5ad",
        out_prefix="batch_check_adult",
        label="Adult",
    )
