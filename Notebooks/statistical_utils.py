"""Statistical helpers for donor-level communication analyses."""

from __future__ import annotations

import hashlib

import numpy as np

from analysis_config import BOOTSTRAP_REPLICATES, RANDOM_SEED


def _interaction_seed(interaction_id: str) -> int:
    digest = hashlib.sha256(interaction_id.encode("utf-8")).digest()
    return (int.from_bytes(digest[:4], "little") + RANDOM_SEED) % (2**32)


def pairwise_effect_sizes(values_a, values_b, interaction_id: str) -> dict:
    """Return interpretable effect sizes and a bootstrap median-difference CI."""
    a = np.asarray(values_a, dtype=float)
    b = np.asarray(values_b, dtype=float)
    if a.ndim != 1 or b.ndim != 1 or not len(a) or not len(b):
        raise ValueError("Both groups must contain one-dimensional observations.")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Effect-size inputs must be finite.")

    rng = np.random.default_rng(_interaction_seed(interaction_id))
    a_boot = a[rng.integers(0, len(a), size=(BOOTSTRAP_REPLICATES, len(a)))]
    b_boot = b[rng.integers(0, len(b), size=(BOOTSTRAP_REPLICATES, len(b)))]
    differences = np.median(a_boot, axis=1) - np.median(b_boot, axis=1)
    ci_low, ci_high = np.quantile(differences, [0.025, 0.975])
    return {
        "median_difference_a_minus_b": float(np.median(a) - np.median(b)),
        "median_difference_ci95_low": float(ci_low),
        "median_difference_ci95_high": float(ci_high),
    }


def rank_biserial_from_u(u_statistic: float, n_a: int, n_b: int) -> float:
    """Rank-biserial correlation oriented as group A minus group B."""
    return float(2 * u_statistic / (n_a * n_b) - 1)


def approximate_spearman_ci(rho: float, n: int) -> tuple[float, float]:
    """Approximate 95% CI using Fisher's z transform; intended for reporting."""
    if n <= 3 or not np.isfinite(rho):
        return (np.nan, np.nan)
    clipped = np.clip(rho, -0.999999, 0.999999)
    z = np.arctanh(clipped)
    half_width = 1.96 / np.sqrt(n - 3)
    return tuple(map(float, np.tanh([z - half_width, z + half_width])))
