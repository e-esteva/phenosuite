"""
circuit_engine.py — Core n-way spatial co-localization engine.

Single source of truth for circuit enrichment analysis, used by:
  - spatial-dynamics CLI (via run-spatial_circuit-enrichment.py)
  - circuit_enrichment Shiny app (via reticulate)

All functions accept and return numpy arrays / plain dicts so they
serialize cleanly across the reticulate bridge.
"""

import numpy as np
from scipy.spatial import cKDTree


def compute_neighborhood_composition(xy, celltypes, radius, circuit_types):
    """
    For every cell, count how many of each circuit member fall within `radius`.

    Parameters
    ----------
    xy : ndarray, shape (n, 2)
        Cell centroid coordinates.
    celltypes : ndarray of str, shape (n,)
        Cell-type label per cell.
    radius : float
        Neighborhood radius in coordinate units.
    circuit_types : list of str
        Ordered list of circuit member cell types.

    Returns
    -------
    comp : ndarray, shape (n, len(circuit_types))
        Count of each circuit member within radius of each cell.
    """
    tree = cKDTree(xy)
    neighbors_list = tree.query_ball_tree(tree, r=radius)

    n = len(xy)
    n_circuit = len(circuit_types)
    ct_to_idx = {ct: i for i, ct in enumerate(circuit_types)}
    comp = np.zeros((n, n_circuit), dtype=np.int32)

    for i in range(n):
        for j in neighbors_list[i]:
            ct = celltypes[j]
            if ct in ct_to_idx:
                comp[i, ct_to_idx[ct]] += 1

    return comp


def circuit_score(comp, method="min_fraction"):
    """
    Score each neighborhood for circuit completeness.

    Parameters
    ----------
    comp : ndarray, shape (n, k)
        Neighborhood composition matrix from compute_neighborhood_composition.
    method : str
        'min_fraction' — bottlenecked by rarest member (strict n-simplex test).
        'geometric_mean' — geometric mean of per-member fractions (softer).

    Returns
    -------
    scores : ndarray, shape (n,)
        Circuit completeness score per cell neighborhood.
    """
    row_totals = np.maximum(comp.sum(axis=1), 1).astype(np.float64)
    fracs = comp / row_totals[:, np.newaxis]

    if method == "min_fraction":
        return fracs.min(axis=1)
    elif method == "geometric_mean":
        log_fracs = np.log(np.maximum(fracs, 1e-10))
        return np.exp(log_fracs.mean(axis=1))
    else:
        raise ValueError(f"Unknown method: {method}")


def circuit_zscore(scores, comp, n_perm=500, seed=42):
    """
    Permutation-based z-score for circuit enrichment.

    Shuffles cell-type labels (by permuting rows of the composition matrix),
    recomputes scores under the null, and compares observed mean to null.

    Parameters
    ----------
    scores : ndarray, shape (n,)
        Observed circuit scores.
    comp : ndarray, shape (n, k)
        Observed composition matrix.
    n_perm : int
        Number of permutations.
    seed : int
        RNG seed for reproducibility.

    Returns
    -------
    result : dict
        z, obs_mean, null_mean, null_sd, p_value
    """
    rng = np.random.default_rng(seed)
    obs_mean = float(scores.mean())
    n = comp.shape[0]

    null_means = np.empty(n_perm)
    for p in range(n_perm):
        perm_comp = comp[rng.permutation(n)]
        perm_totals = np.maximum(perm_comp.sum(axis=1), 1).astype(np.float64)
        perm_fracs = perm_comp / perm_totals[:, np.newaxis]
        perm_scores = perm_fracs.min(axis=1)
        null_means[p] = perm_scores.mean()

    null_mean = float(null_means.mean())
    null_sd = float(null_means.std())
    z = (obs_mean - null_mean) / max(null_sd, 1e-10)
    p_value = float((null_means >= obs_mean).mean())

    return {
        "z": z,
        "obs_mean": obs_mean,
        "null_mean": null_mean,
        "null_sd": null_sd,
        "p_value": p_value,
    }


def threshold_sweep(scores, thresholds):
    """
    For each threshold, compute n_positive and frac_positive.

    Parameters
    ----------
    scores : ndarray, shape (n,)
    thresholds : list or ndarray of float

    Returns
    -------
    results : list of dict
        Each dict has: threshold, n_positive, frac_positive, mean_score
    """
    thresholds = np.asarray(thresholds)
    results = []
    n = len(scores)
    for t in thresholds:
        pos = scores >= t
        n_pos = int(pos.sum())
        results.append({
            "threshold": float(t),
            "n_positive": n_pos,
            "frac_positive": n_pos / max(n, 1),
            "mean_score": float(scores[pos].mean()) if n_pos > 0 else 0.0,
        })
    return results
