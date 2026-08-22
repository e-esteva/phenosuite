"""
circuit_engine.py — Core n-way spatial co-localization engine.

Hybrid numpy / scipy.sparse / JAX implementation.

Dispatch strategy:
  always               → scipy.sparse matmul for neighborhood composition
  n < LARGE_N (100 k)  → numpy for permutation test
  n ≥ LARGE_N          → JAX vmap + JIT on CPU via XLA

Used by:
  - spatial-dynamics CLI (via run-spatial_circuit-enrichment.py)
  - circuit_enrichment Shiny app (via reticulate)

All functions accept and return numpy arrays / plain dicts so they
serialize cleanly across the reticulate bridge.
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix

LARGE_N = 100_000
_PERM_CHUNK = 20   # permutations per vmap batch; keeps peak memory ≈ chunk × n × k × 4 B

_jax_cache = {}
_jit_cache = {}


def _get_jax():
    if "jax" not in _jax_cache:
        import jax
        import jax.numpy as jnp
        jax.config.update("jax_platform_name", "cpu")
        _jax_cache["jax"] = jax
        _jax_cache["jnp"] = jnp
    return _jax_cache["jax"], _jax_cache["jnp"]


def _get_perm_fn(method):
    """
    Return the cached jit(vmap(single_permutation)) function for `method`.

    Null model: independently shuffles each circuit member's neighbor-count
    column across cells. That breaks the cross-member co-occurrence this test
    is meant to detect while preserving each member's own marginal
    distribution of local abundance — permuting whole rows (the old approach)
    doesn't touch a row-order-invariant statistic like this one at all, so it
    can't produce a real null distribution.
    """
    cache_key = f"perm_{method}"
    if cache_key not in _jit_cache:
        jax, jnp = _get_jax()

        def _one_perm(key, comp_j):
            n, k = comp_j.shape
            col_keys = jax.random.split(key, k)
            cols = [jax.random.permutation(col_keys[j], comp_j[:, j]) for j in range(k)]
            pc = jnp.stack(cols, axis=1)
            totals = jnp.maximum(pc.sum(axis=1, keepdims=True), 1.0)
            fracs = pc / totals
            if method == "min_fraction":
                return fracs.min(axis=1).mean()
            elif method == "geometric_mean":
                log_fracs = jnp.log(jnp.maximum(fracs, 1e-10))
                return jnp.exp(log_fracs.mean(axis=1)).mean()
            else:
                raise ValueError(f"Unknown method: {method}")

        _jit_cache[cache_key] = jax.jit(jax.vmap(_one_perm, in_axes=(0, None)))
    return _jit_cache[cache_key]


# ============================================================================
# PUBLIC API
# ============================================================================

def compute_neighborhood_composition(xy, celltypes, radius, circuit_types):
    """
    For every cell, count how many of each circuit member fall within `radius`.

    Builds a sparse adjacency matrix in vectorized COO format, then multiplies
    by a one-hot cell-type indicator to produce the composition counts in a
    single scipy.sparse @ dense call — no Python loops over cells.

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
    comp : ndarray, shape (n, len(circuit_types)), int32
        Count of each circuit member within radius of each cell.
    """
    tree = cKDTree(xy)
    neighbors_list = tree.query_ball_tree(tree, r=radius)

    n = len(xy)
    n_circuit = len(circuit_types)
    ct_to_idx = {ct: i for i, ct in enumerate(circuit_types)}

    # Vectorized COO construction — avoids inner Python loop
    lengths = np.fromiter((len(nb) for nb in neighbors_list), dtype=np.int32, count=n)
    row_idx = np.repeat(np.arange(n, dtype=np.int32), lengths)
    if lengths.sum() > 0:
        col_idx = np.concatenate(neighbors_list).astype(np.int32)
    else:
        col_idx = np.empty(0, dtype=np.int32)

    A = csr_matrix(
        (np.ones(len(row_idx), dtype=np.float32), (row_idx, col_idx)),
        shape=(n, n),
    )

    # One-hot cell-type indicator, shape (n, n_circuit)
    ct_indices = np.array([ct_to_idx.get(ct, -1) for ct in celltypes], dtype=np.int32)
    valid = np.where(ct_indices >= 0)[0]
    B = np.zeros((n, n_circuit), dtype=np.float32)
    B[valid, ct_indices[valid]] = 1.0

    # Single sparse @ dense matmul replaces the double Python loop
    return np.asarray(A @ B, dtype=np.int32)


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


def _score_from_comp(pc, method):
    """Mean circuit-completeness score over all rows of a composition matrix — the same math as circuit_score(), reused so the null distribution is scored identically to the observed data."""
    totals = np.maximum(pc.sum(axis=1), 1).astype(np.float64)
    fracs = pc / totals[:, np.newaxis]
    if method == "min_fraction":
        return fracs.min(axis=1).mean()
    elif method == "geometric_mean":
        log_fracs = np.log(np.maximum(fracs, 1e-10))
        return np.exp(log_fracs.mean(axis=1)).mean()
    else:
        raise ValueError(f"Unknown method: {method}")


def circuit_zscore(scores, comp, n_perm=500, seed=42, method="min_fraction"):
    """
    Permutation-based z-score for circuit enrichment.

    Null model: independently permutes each circuit member's neighbor-count
    column across cells (see _get_perm_fn / _zscore_numpy), then rescores
    with the same `method` used for the observed `scores` so obs_mean and
    null_mean are directly comparable.

    Routes to JAX vmap path (XLA-threaded, chunked) for n >= LARGE_N,
    and to a plain numpy loop otherwise.

    Parameters
    ----------
    scores : ndarray, shape (n,)
        Observed circuit scores (from circuit_score(comp, method)).
    comp : ndarray, shape (n, k)
        Observed composition matrix.
    n_perm : int
        Number of permutations.
    seed : int
        RNG seed for reproducibility.
    method : str
        Scoring method used for `scores` — 'min_fraction' or 'geometric_mean' —
        so the null distribution is scored the same way.

    Returns
    -------
    result : dict
        z, obs_mean, null_mean, null_sd, p_value
    """
    obs_mean = float(scores.mean())
    n = comp.shape[0]

    null_means = (_zscore_jax(comp, n_perm, seed, method) if n >= LARGE_N
                  else _zscore_numpy(comp, n_perm, seed, method))

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

    Vectorized over all thresholds simultaneously.

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
    n = len(scores)
    # Broadcast comparison: (n, 1) >= (1, T) → (n, T)
    positive = scores[:, np.newaxis] >= thresholds[np.newaxis, :]
    n_pos = positive.sum(axis=0)

    results = []
    for t_idx, t in enumerate(thresholds):
        k = int(n_pos[t_idx])
        results.append({
            "threshold": float(t),
            "n_positive": k,
            "frac_positive": k / max(n, 1),
            "mean_score": float(scores[positive[:, t_idx]].mean()) if k > 0 else 0.0,
        })
    return results


# ============================================================================
# INTERNAL HELPERS
# ============================================================================

def _zscore_numpy(comp, n_perm, seed, method):
    """
    Independently permutes each circuit member's column across cells (rather
    than whole rows) so the co-occurrence structure the test targets is
    actually destroyed under the null, while each member's own marginal
    abundance distribution is preserved.
    """
    rng = np.random.default_rng(seed)
    n, k = comp.shape
    null_means = np.empty(n_perm)
    for p in range(n_perm):
        pc = np.column_stack([rng.permutation(comp[:, j]) for j in range(k)])
        null_means[p] = _score_from_comp(pc, method)
    return null_means


def _zscore_jax(comp, n_perm, seed, method):
    """
    Chunked vmap over permutations.

    Processes _PERM_CHUNK permutations per XLA call so peak memory stays
    bounded at _PERM_CHUNK × n × k × 4 bytes regardless of n_perm.
    The jit(vmap(fn)) kernel is cached per-method across calls via _jit_cache.
    """
    jax, jnp = _get_jax()
    perm_fn = _get_perm_fn(method)
    comp_j = jnp.array(comp, dtype=jnp.float32)

    all_keys = jax.random.split(jax.random.PRNGKey(seed), n_perm)
    chunks = []
    for start in range(0, n_perm, _PERM_CHUNK):
        chunk_keys = all_keys[start : start + _PERM_CHUNK]
        chunks.append(np.array(perm_fn(chunk_keys, comp_j)))

    return np.concatenate(chunks)
