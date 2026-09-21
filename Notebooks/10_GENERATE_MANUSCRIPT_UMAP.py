"""Generate the combined pediatric/adult cell-type UMAP for the manuscript."""

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from publication_utils import save_publication_figure


BASE_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = BASE_DIR / "Results" / "manuscript_figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH = OUTPUT_DIR / "figure1_pediatric_adult_celltype_umap.png"

INPUTS = {
    "A. Pediatric Cohort": (
        BASE_DIR / "child_individual_h5ad" / "adata_child_final_liana.h5ad"
    ),
    "B. Adult Cohort": (
        BASE_DIR / "adult_individual_h5ad" / "adata_adult_final_liana.h5ad"
    ),
}

CELL_TYPE_ORDER = [
    "B Cells",
    "CD4 T Cells",
    "CD8 T Cells",
    "NK Cells",
    "Classical Monocytes",
    "Non-classical Monocytes",
    "Plasma Cells",
    "pDCs",
    "cDCs",
    "IFN-stimulated T Cells",
]

# The same cell type always receives the same color in both panels.
COLORS = {
    "B Cells": "#1f77b4",
    "CD4 T Cells": "#ff7f0e",
    "CD8 T Cells": "#2ca02c",
    "NK Cells": "#9467bd",
    "Classical Monocytes": "#d62728",
    "Non-classical Monocytes": "#8c564b",
    "Plasma Cells": "#e377c2",
    "pDCs": "#f2a541",
    "cDCs": "#17becf",
    "IFN-stimulated T Cells": "#7f7f7f",
}


def load_umap(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read only the final UMAP coordinates and cell-type labels."""
    adata = ad.read_h5ad(path, backed="r")
    try:
        coordinates = np.asarray(adata.obsm["X_umap"]).copy()
        labels = adata.obs["cell_type"].astype(str).to_numpy(copy=True)
    finally:
        adata.file.close()
    return coordinates, labels


def main() -> None:
    datasets = {title: load_umap(path) for title, path in INPUTS.items()}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2))

    for ax, (title, (coordinates, labels)) in zip(axes, datasets.items()):
        # Plot larger populations first and rare populations last so rare cells remain visible.
        counts = {cell_type: int(np.sum(labels == cell_type)) for cell_type in CELL_TYPE_ORDER}
        draw_order = sorted(CELL_TYPE_ORDER, key=lambda cell_type: counts[cell_type], reverse=True)
        for cell_type in draw_order:
            selected = labels == cell_type
            if not selected.any():
                continue
            ax.scatter(
                coordinates[selected, 0],
                coordinates[selected, 1],
                s=0.45,
                c=COLORS[cell_type],
                alpha=0.72,
                linewidths=0,
                rasterized=True,
            )
        ax.set_title(title, fontsize=17, fontweight="bold", loc="left")
        ax.set_xlabel("UMAP 1", fontsize=15)
        ax.set_ylabel("UMAP 2", fontsize=15)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=8,
            markerfacecolor=COLORS[cell_type],
            markeredgecolor="none",
            label=cell_type,
        )
        for cell_type in CELL_TYPE_ORDER
    ]
    fig.legend(
        handles=legend_handles,
        title="Cell Type",
        loc="center left",
        bbox_to_anchor=(0.82, 0.5),
        frameon=False,
        fontsize=13,
        title_fontsize=14,
    )
    fig.subplots_adjust(left=0.05, right=0.81, bottom=0.10, top=0.92, wspace=0.12)
    save_publication_figure(OUTPUT_PATH, fig)
    plt.close(fig)
    print(f"Wrote 600-DPI manuscript UMAP to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
