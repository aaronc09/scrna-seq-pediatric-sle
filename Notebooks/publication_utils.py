"""Shared reproducibility and publication-output helpers.

All notebooks import this module.  It deliberately avoids project-specific biology;
its only responsibilities are runtime validation and safe, deterministic output.
"""

from __future__ import annotations

import platform
import shutil
import json
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt


REQUIRED_PYTHON = (3, 12, 10)
FIGURE_DPI = 600
LEGACY_FIGURES_BY_NOTEBOOK = {
    "01B_GSE135779_CHILD_CELLTYPE_LABELLING": (
        "cluster3_child_ribo_cd8_umap.png",
        "cluster3_child_subcluster_umap.png",
    ),
    "01E_GSE135779_ADULT_CELLTYPE_LABELLING": (
        "FINAL_umap_child_adult_celltype.png",
        "FINAL_umap_child_adult_celltype.pdf",
    ),
}
LEGACY_OUTPUT_DIRS_BY_NOTEBOOK = {
    "01A_GSE135779_CHILD_PREPROCESSING": "1A_GSE135779_CHILD_PREPROCESSING",
    "01B_GSE135779_CHILD_CELLTYPE_LABELLING": "1B_GSE135779_CHILD_CELLTYPE_LABELLING",
    "01D_GSE135779_ADULT_PREPROCESSING": "1D_GSE135779_ADULT_PREPROCESSING",
    "01E_GSE135779_ADULT_CELLTYPE_LABELLING": "1E_GSE135779_ADULT_CELLTYPE_LABELLING",
}
LEGACY_INTERMEDIATES_BY_NOTEBOOK = {
    "01B_GSE135779_CHILD_CELLTYPE_LABELLING": (
        "child_individual_h5ad/adata_child_cluster3_subclustered.h5ad",
    ),
}


def require_python_version() -> None:
    """Fail early when the analysis is run with a different Python version."""
    running = tuple(map(int, platform.python_version_tuple()))
    if running != REQUIRED_PYTHON:
        required = ".".join(map(str, REQUIRED_PYTHON))
        raise RuntimeError(
            f"This analysis requires Python {required}; found {platform.python_version()}."
        )


def configure_publication_notebook(base_dir: str | Path, notebook_stem: str) -> Path:
    """Configure consistent plotting and return this notebook's figure directory."""
    require_python_version()
    figure_dir = Path(base_dir) / "Notebook_Outputs" / notebook_stem
    figure_dir.mkdir(parents=True, exist_ok=True)
    # This directory is owned by one notebook, so removing its old figures is safe.
    # Tables are deliberately preserved; analytical result files are overwritten by
    # their producing cells and may be required as inputs by later notebooks.
    for stale_figure in figure_dir.glob("*.png"):
        stale_figure.unlink()
    output_root = Path(base_dir) / "Notebook_Outputs"
    legacy_directory_name = LEGACY_OUTPUT_DIRS_BY_NOTEBOOK.get(notebook_stem)
    if legacy_directory_name:
        legacy_directory = output_root / legacy_directory_name
        if legacy_directory.is_dir():
            shutil.rmtree(legacy_directory)
    for filename in LEGACY_FIGURES_BY_NOTEBOOK.get(notebook_stem, ()):
        legacy_path = output_root / filename
        if legacy_path.is_file():
            legacy_path.unlink()
    for relative_path in LEGACY_INTERMEDIATES_BY_NOTEBOOK.get(notebook_stem, ()):
        legacy_intermediate = Path(base_dir) / relative_path
        if legacy_intermediate.is_file():
            legacy_intermediate.unlink()
    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": FIGURE_DPI,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
            "font.family": "sans-serif",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    # Scanpy uses its own figure directory and dpi_save setting.
    try:
        import scanpy as sc

        sc.settings.figdir = figure_dir
        sc.settings.set_figure_params(dpi_save=FIGURE_DPI, facecolor="white")
    except ImportError:
        # Some metadata-only notebooks do not otherwise need Scanpy.
        pass
    provenance_dir = Path(base_dir) / "Results" / "provenance"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    packages = {}
    for package in ("anndata", "harmonypy", "liana", "numpy", "pandas", "scanpy",
                    "scipy", "scrublet", "statsmodels"):
        try:
            packages[package] = version(package)
        except PackageNotFoundError:
            packages[package] = None
    manifest = {
        "notebook": notebook_stem,
        "run_started_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": packages,
    }
    (provenance_dir / f"{notebook_stem}.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return figure_dir


def remove_owned_outputs(paths: Iterable[str | Path]) -> None:
    """Remove only explicitly declared files owned by the current notebook."""
    for raw_path in paths:
        path = Path(raw_path)
        if path.exists():
            if not path.is_file():
                raise IsADirectoryError(
                    f"Refusing to remove a directory as an owned output: {path}"
                )
            path.unlink()


def remove_owned_intermediates(base_dir: str | Path, relative_paths: Iterable[str]) -> None:
    """Remove declared disposable checkpoints without permitting path traversal."""
    base = Path(base_dir).resolve()
    for relative_path in relative_paths:
        candidate = (base / relative_path).resolve()
        if candidate == base or base not in candidate.parents:
            raise ValueError(f"Refusing to remove path outside the project: {candidate}")
        if candidate.is_dir():
            shutil.rmtree(candidate)
        elif candidate.is_file():
            candidate.unlink()


def save_publication_figure(path: str | Path, fig=None) -> Path:
    """Save a figure as a 600-DPI PNG, replacing an older copy atomically."""
    output_path = Path(path)
    if output_path.suffix.lower() != ".png":
        raise ValueError(f"Publication figures must use PNG format: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(f".{output_path.stem}.tmp.png")
    target_figure = fig if fig is not None else plt.gcf()
    target_figure.savefig(
        temporary_path,
        dpi=FIGURE_DPI,
        bbox_inches="tight",
        facecolor="white",
    )
    temporary_path.replace(output_path)
    return output_path
