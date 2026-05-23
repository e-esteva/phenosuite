"""
pwlo-es-pt.py — Pairwise log-odds interaction engine.

Uses scipy cKDTree for O(N log M) neighbor counting with negligible memory
overhead. query_ball_point with return_length=True returns only counts per
point, never materialising index lists or distance matrices.

Effect size (optional): KS statistic comparing actual vs random neighbor
distances, subsampled to keep memory bounded.
"""


def pairwise_logOdds(spatial_obj, out_dir, label, draw=False, resolution=0.3774, p1=3, p2=30, compute_effect_size=False):
    import numpy as np
    import pandas as pd
    from scipy.stats import ks_2samp
    from scipy.spatial import cKDTree
    import os

    print(spatial_obj)
    print(compute_effect_size)

    set_sizes  = spatial_obj['cluster'].value_counts()
    celltypes  = set_sizes.keys()

    p2         = p2 + p1
    p1_scaled  = p1 / resolution
    p2_scaled  = p2 / resolution

    # Number of ref cells sampled per pair for the KS distance comparison.
    # Kept small so the distance sample vectors are tractable without cdist.
    KS_REF_SAMPLE  = 2000
    KS_DIST_SAMPLE = 10000

    global_logodds      = []
    global_probabilities = []
    global_effect_sizes  = []

    for i in range(len(celltypes)):

        ref_      = celltypes[i]
        ref_data  = spatial_obj[['x', 'y']][spatial_obj['cluster'] == ref_].values

        neighbors_i = []
        counts_i    = []
        blocks_ks   = []

        for j in range(len(celltypes)):

            query_     = celltypes[j]
            query_data = spatial_obj[['x', 'y']][spatial_obj['cluster'] == query_].values

            print(f'{ref_} | N= {set_sizes[ref_]} : {query_} | N= {set_sizes[query_]}')

            # ------------------------------------------------------------------
            # Core neighbor count — O(N log M) time, O(N + M) memory.
            # Build a KD-tree on query cells, then count how many fall inside
            # the annulus [p1_scaled, p2_scaled] around each ref cell.
            # return_length=True means only the count is returned per point,
            # never the full list of indices, keeping peak RAM negligible.
            # ------------------------------------------------------------------
            query_tree  = cKDTree(query_data)

            outer = query_tree.query_ball_point(ref_data, r=p2_scaled, return_length=True, workers=-1)
            inner = query_tree.query_ball_point(ref_data, r=p1_scaled, return_length=True, workers=-1)

            annulus = outer - inner                          # per-ref-cell neighbor count

            neighbors_i.append(int(np.sum(annulus)))
            counts_i.append(int(np.sum(annulus > 0)))

            # ------------------------------------------------------------------
            # Optional KS effect size.
            # Instead of materialising an N×M distance matrix we:
            #   1. Subsample KS_REF_SAMPLE ref cells.
            #   2. Retrieve actual neighbor distances (within p2_scaled) via the
            #      tree — O(ref_sample × mean_neighbors) rather than O(N×M).
            #   3. Repeat against a size-matched random draw from the full data.
            #   4. Downsample both distance vectors to KS_DIST_SAMPLE before
            #      the KS call so the test stays fast regardless of density.
            # ------------------------------------------------------------------
            if compute_effect_size:

                n_ref_sample = min(len(ref_data), KS_REF_SAMPLE)
                sample_idx   = np.random.choice(len(ref_data), size=n_ref_sample, replace=False)
                ref_sample   = ref_data[sample_idx]

                # Actual distances to query neighbours within the search radius
                neighbor_idx_lists = query_tree.query_ball_point(ref_sample, r=p2_scaled, workers=-1)
                actual_dists = np.concatenate([
                    np.linalg.norm(query_data[idx] - ref_sample[k], axis=1)
                    for k, idx in enumerate(neighbor_idx_lists)
                    if len(idx) > 0
                ]) if any(len(idx) > 0 for idx in neighbor_idx_lists) else np.array([])

                # Random-baseline distances: same query size, random positions
                rand_data   = spatial_obj[['x', 'y']].sample(len(query_data)).values
                rand_tree   = cKDTree(rand_data)
                rand_idx_lists = rand_tree.query_ball_point(ref_sample, r=p2_scaled, workers=-1)
                random_dists = np.concatenate([
                    np.linalg.norm(rand_data[idx] - ref_sample[k], axis=1)
                    for k, idx in enumerate(rand_idx_lists)
                    if len(idx) > 0
                ]) if any(len(idx) > 0 for idx in rand_idx_lists) else np.array([])

                if len(actual_dists) > 1 and len(random_dists) > 1:
                    n_ks    = min(len(actual_dists), len(random_dists), KS_DIST_SAMPLE)
                    ad_s    = np.random.choice(actual_dists,  size=n_ks, replace=False)
                    rd_s    = np.random.choice(random_dists,  size=n_ks, replace=False)
                    ks_stat = ks_2samp(ad_s, rd_s, alternative='greater').statistic
                else:
                    ks_stat = 0.0

                blocks_ks.append(ks_stat)
                print(f'effect size: {ks_stat}')

        # ----------------------------------------------------------------------
        # Log-odds row for this ref type.
        # Guard against zero total neighbors (all-zero row → log(0)) with a
        # small epsilon so we produce -inf gracefully rather than crashing.
        # ----------------------------------------------------------------------
        total_neighbors = np.sum(np.array(neighbors_i))

        log_odds_i = np.array([
            np.log(
                (neighbors_i[x] / (total_neighbors + 1e-12)) /
                (set_sizes[celltypes[x]] / spatial_obj.shape[0])
            )
            for x in range(len(neighbors_i))
        ])
        print('log_odds_i: ' + str(log_odds_i))
        global_logodds.append(log_odds_i)

        global_probabilities.append(np.array(counts_i) / set_sizes[ref_])

        if compute_effect_size:
            global_effect_sizes.append(blocks_ks)

    # --------------------------------------------------------------------------
    # Assemble and write output matrices
    # --------------------------------------------------------------------------
    global_logodds_df           = pd.DataFrame(global_logodds).T
    global_logodds_df.columns   = celltypes
    global_logodds_df.index     = celltypes
    global_logodds_df.to_csv(str(out_dir) + '/' + str(label) + '-logOdds_matrix.csv')

    global_probabilities_df         = pd.DataFrame(global_probabilities).T
    global_probabilities_df.columns = celltypes
    global_probabilities_df.index   = celltypes
    global_probabilities_df.to_csv(str(out_dir) + '/' + str(label) + '-probabilities_matrix.csv')

    if compute_effect_size:
        global_effect_sizes_df         = pd.DataFrame(global_effect_sizes).T
        global_effect_sizes_df.columns = celltypes
        global_effect_sizes_df.index   = celltypes
        global_effect_sizes_df.to_csv(str(out_dir) + '/' + str(label) + '-KS-effect_sizes_matrix.csv')

    return global_logodds_df
