# NeighborhoodR ── Spatial Neighborhood Analysis (memory-optimised)
# PhenoSuite | neighborhood_analysis/app.R
#
# Python backend (reticulate + sklearn/scipy/numpy) is used for:
#   • KNN  — sklearn NearestNeighbors (parallel, ball_tree)
#   • Niche matrix — scipy.sparse scatter-add (float32, ~50 % RAM vs R float64)
#   • Clustering  — sklearn MiniBatchKMeans (memory-efficient for M-scale cell counts)
#   • LOO sweep   — fully in Python; only small summary vectors returned to R
# R pure-R fallback (RANN + kmeans) activates automatically when sklearn is absent.
# ─────────────────────────────────────────────────────────────────────────────

# Auto-install any packages missing from the local environment (no-op on the
# server where everything is pre-baked into the Docker image).
.ensure <- function(...) {
  pkgs <- c(...)
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}
.ensure("shiny", "ggplot2", "dplyr", "tidyr", "purrr",
        "DT", "scales", "shinycssloaders", "glue", "jsonlite", "RANN")

suppressPackageStartupMessages({
  library(shiny)
  library(reticulate)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(RANN)          # R fallback KNN
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(DT)
  library(scales)
  library(shinycssloaders)
  library(glue)
  library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b

prov_path <- '/srv/shiny-server/phenomenalist/utils/provenance.R'
if (file.exists(prov_path)) source(prov_path)

options(shiny.maxRequestSize = 4 * 500 * 1024^2)

# ─── Python setup ────────────────────────────────────────────────────────────
# All compute-heavy work runs in Python if sklearn is available.
# The niche matrix is kept as a Python float32 object through the full pipeline.

PY_CODE <- r"(
import numpy as np
from scipy.sparse import coo_matrix as _coo
from sklearn.neighbors import NearestNeighbors as _NN
from sklearn.cluster import MiniBatchKMeans as _MBK
import gc as _gc, os as _os

_PSUTIL = False
try:
    import psutil as _ps
    _PSUTIL = True
except ImportError:
    pass

def memory_mb():
    if _PSUTIL:
        return round(_ps.Process(_os.getpid()).memory_info().rss / 1024**2, 1)
    return None

def deep_clean():
    """
    Full heap release: Python GC + ask glibc to return freed pages to the OS.
    malloc_trim(0) is a no-op on non-Linux but harmless.
    """
    _gc.collect()
    try:
        import ctypes
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass
    return True

def build_pooled_niche_matrix(coords_list, ct_enc_list, n_cts, k1):
    """
    Build (total_cells x n_cts) float32 niche matrix across all samples.
    coords_list  : list of (n_i x 2) float matrices (R list of matrices)
    ct_enc_list  : list of (n_i,) 0-indexed int vectors   (R list of integers)
    n_cts        : int, global celltype count
    k1           : int, KNN neighbours

    Keeps everything in float32 and processes samples one at a time to
    minimise peak memory.
    """
    n_cts = int(n_cts)
    k1    = int(k1)

    # first pass: compute total cells so we can pre-allocate
    sizes = [np.asarray(c, dtype=np.float32).shape[0] for c in coords_list]
    total = int(sum(sizes))
    niche = np.zeros((total, n_cts), dtype=np.float32)

    ptr = 0
    for coords_r, ct_enc_r in zip(coords_list, ct_enc_list):
        coords = np.asarray(coords_r, dtype=np.float32)
        ct_enc = np.asarray(ct_enc_r, dtype=np.int32).ravel()
        n      = coords.shape[0]
        k      = min(k1, n - 1)

        nbrs = _NN(n_neighbors=k + 1, algorithm='ball_tree',
                   n_jobs=-1, metric='euclidean')
        nbrs.fit(coords)
        _, idx = nbrs.kneighbors(coords)
        idx = idx[:, 1:]          # (n x k) – drop self

        # sparse scatter-add: much faster than np.add.at for large n
        row  = np.repeat(np.arange(n, dtype=np.int32), k)
        col  = ct_enc[idx.ravel()]
        vals = np.ones(len(row), dtype=np.float32)
        block = _coo((vals, (row, col)), shape=(n, n_cts)).toarray()
        block /= k

        niche[ptr:ptr + n] = block
        ptr += n

        del coords, ct_enc, idx, row, col, vals, block
        _gc.collect()

    return niche   # float32 — stays in Python heap

def loo_stability_sweep(niche_mat, ct_enc_r, samp_enc_r, n_samples,
                        k_sweep_r, loo_n, loo_mode, agg, n_cts,
                        sample_group_r=None):
    """
    Full LOO stability sweep in Python.
    niche_mat : float32 numpy array (already in Python memory)
    ct_enc_r  : 0-indexed celltype vector (R integer vector)
    samp_enc_r: 0-indexed sample vector   (R integer vector)
    sample_group_r : 0-indexed group-per-sample vector (R integer vector),
        required when loo_mode == 'group' — e.g. timepoint, from the mapped
        condition. Holds out loo_n samples from *every* group each fold,
        instead of drawing loo_n/pct from the pooled sample list, so a
        stratified structure (e.g. N timepoints x R replicates) actually
        gets tested as intended rather than approximated by picking a total
        count and hoping the random draw lands evenly across groups.
    Returns dict with 'k' and 'stability' lists.
    """
    nm       = niche_mat                        # reference, no copy
    ct_enc   = np.asarray(ct_enc_r,  dtype=np.int32).ravel()
    samp_enc = np.asarray(samp_enc_r, dtype=np.int32).ravel()
    k_sweep  = [int(k) for k in k_sweep_r]
    n_s      = int(n_samples)
    n_cts    = int(n_cts)
    loo_n    = float(loo_n)

    rng = np.random.default_rng(0)
    n_iter = min(n_s, 20)

    if loo_mode == 'group':
        sample_group = np.asarray(sample_group_r, dtype=np.int32).ravel()
        groups       = np.unique(sample_group)
        group_idx    = {g: np.where(sample_group == g)[0] for g in groups}
        min_group_sz = min(len(idx) for idx in group_idx.values())
        # Need >=1 training sample left in the smallest group.
        n_hold_per_group = max(1, min(int(loo_n), min_group_sz - 1))
        if n_hold_per_group < 1 or len(groups) < 2:
            return {'k': k_sweep, 'stability': [0.0] * len(k_sweep)}
        loo_sets = [
            np.concatenate([rng.choice(idx, n_hold_per_group, replace=False)
                             for idx in group_idx.values()])
            for _ in range(n_iter)
        ]
    else:
        n_hold = max(1, min(int(loo_n), n_s - 1)) if loo_mode == 'count' \
                 else max(1, int(np.floor(n_s * loo_n / 100.0)))

        # Can't hold out if we'd remove all samples — return zero instability
        if n_hold >= n_s:
            return {'k': k_sweep, 'stability': [0.0] * len(k_sweep)}

        loo_sets = [rng.choice(n_s, n_hold, replace=False) for _ in range(n_iter)]

    agg_fn = np.median if agg == 'median' else np.mean

    def _nh_freq(assign, ct, k2):
        freq = np.zeros((k2, n_cts), dtype=np.float32)
        for j in range(k2):
            m = assign == j
            if not m.any():
                continue
            bc = np.bincount(ct[m], minlength=n_cts).astype(np.float32)
            freq[j] = bc / m.sum()
        return freq

    def _nearest(data, centers):
        chunk = 8192
        out   = np.empty(data.shape[0], dtype=np.int32)
        for s in range(0, data.shape[0], chunk):
            e    = min(s + chunk, data.shape[0])
            d    = data[s:e, np.newaxis, :] - centers[np.newaxis, :, :]
            out[s:e] = np.argmin(np.einsum('ijk,ijk->ij', d, d), axis=1)
        return out

    scores = []
    for k2 in k_sweep:
        deltas = []
        for held in loo_sets:
            held_set   = set(held.tolist())
            train_mask = np.array([s not in held_set for s in samp_enc.tolist()])
            test_mask  = ~train_mask
            train_data = nm[train_mask]
            if len(train_data) < k2:
                continue
            km = _MBK(n_clusters=k2, n_init=5, random_state=42,
                      batch_size=min(4096, len(train_data)),
                      max_iter=100, max_no_improvement=10)
            try:
                tr_assign = km.fit_predict(train_data)
            except Exception:
                continue

            all_assign = np.empty(len(nm), dtype=np.int32)
            all_assign[train_mask] = tr_assign
            if test_mask.any():
                all_assign[test_mask] = _nearest(
                    nm[test_mask], km.cluster_centers_.astype(np.float32))

            ft = _nh_freq(tr_assign,   ct_enc[train_mask], k2)
            ff = _nh_freq(all_assign,  ct_enc,             k2)
            deltas.append(float(agg_fn(np.abs(ff - ft))))

        scores.append(float(np.mean(deltas)) if deltas else float('nan'))
    _gc.collect()
    return {'k': k_sweep, 'stability': scores}

def final_kmeans(niche_mat, k2, n_init=10, random_state=42):
    """MiniBatchKMeans final assignment (1-indexed labels)."""
    k2 = int(k2)
    bs = min(4096, max(k2 * 10, 1024))
    km = _MBK(n_clusters=k2, n_init=int(n_init), random_state=int(random_state),
              batch_size=bs, max_iter=300, max_no_improvement=30)
    labels = (km.fit_predict(niche_mat).astype(np.int32) + 1).tolist()
    centers = km.cluster_centers_.tolist()
    _gc.collect()
    return {'labels': labels, 'centers': centers}
)"

PY_AVAILABLE <- tryCatch({
  reticulate::py_module_available("sklearn") &&
  reticulate::py_module_available("scipy")   &&
  reticulate::py_module_available("numpy")
}, error = function(e) FALSE)

if (PY_AVAILABLE) {
  reticulate::py_run_string(PY_CODE)
  message("NeighborhoodR: Python backend active (sklearn ",
          reticulate::py_eval("__import__('sklearn').__version__"), ")")
} else {
  message("NeighborhoodR: Python unavailable — using R backend (RANN + kmeans)")
}

# ─── Memory helpers ───────────────────────────────────────────────────────────

mem_mb_r <- function() {
  tryCatch({
    m <- gc(verbose = FALSE, reset = FALSE)
    # Column 2 is the "used (Mb)" column; avoid selecting by name because gc()
    # has three columns all named "(Mb)" which causes sum() to triple-count.
    round(sum(m[, 2L]), 1)
  }, error = function(e) NA_real_)
}

mem_mb_py <- function() {
  if (!PY_AVAILABLE) return(NULL)
  tryCatch({
    val <- reticulate::py$memory_mb()
    # Python returns None when psutil is absent; reticulate maps None -> NULL,
    # and as.numeric(NULL) gives numeric(0) which breaks is.na() checks.
    if (is.null(val) || length(val) == 0L) NULL else as.numeric(val)
  }, error = function(e) NULL)
}

mem_label <- function() {
  r_mb <- mem_mb_r()
  p_mb <- mem_mb_py()
  r_ok  <- isTRUE(!is.na(r_mb))
  py_ok <- isTRUE(length(p_mb) == 1L && !is.na(p_mb))
  if (py_ok && r_ok) {
    sprintf("Process: %.0f MB | R heap: %.0f MB | Python: active", p_mb, r_mb)
  } else if (r_ok) {
    sprintf("R heap: %.0f MB", r_mb)
  } else "—"
}

# ─── Colour palette ──────────────────────────────────────────────────────────
PALETTE_20 <- c(
  "#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
  "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
  "#d4e157","#26c6da","#ab47bc","#7e57c2","#66bb6a",
  "#ffa726","#ef5350","#42a5f5","#26a69a","#c0ca33"
)
pal <- function(n) {
  if (n <= length(PALETTE_20)) PALETTE_20[seq_len(n)]
  else colorRampPalette(PALETTE_20)(n)
}

# ─── R-backend helpers (fallbacks when Python unavailable) ───────────────────

.r_build_niche_matrix <- function(spe, ct_col, k1) {
  coords   <- spatialCoords(spe)
  ctypes   <- as.character(colData(spe)[[ct_col]])
  cts_uniq <- sort(unique(ctypes))
  n_cells  <- ncol(spe)
  k1       <- min(k1, n_cells - 1L)
  nn       <- RANN::nn2(coords, k = k1 + 1L)
  nbr_idx  <- nn$nn.idx[, -1L, drop = FALSE]
  mat <- matrix(0.0, nrow = n_cells, ncol = length(cts_uniq),
                dimnames = list(colnames(spe), cts_uniq))
  for (i in seq_len(n_cells)) {
    tab       <- table(factor(ctypes[nbr_idx[i, ]], levels = cts_uniq))
    mat[i, ]  <- as.numeric(tab) / k1
  }
  mat
}

.r_assign_to_centers <- function(data, centers) {
  apply(data, 1L, function(row)
    which.min(apply(centers, 1L, function(ctr) sum((row - ctr)^2))))
}

compute_nh_freq <- function(assignments, ctypes, k2, all_cts) {
  freq <- matrix(0.0, nrow = k2, ncol = length(all_cts),
                 dimnames = list(paste0("N", seq_len(k2)), all_cts))
  for (j in seq_len(k2)) {
    idx <- which(assignments == j)
    if (!length(idx)) next
    tab        <- table(factor(ctypes[idx], levels = all_cts))
    freq[j, ]  <- as.numeric(tab) / length(idx)
  }
  freq
}

.r_loo_stability_sweep <- function(niche_mat, ctypes_vec, sample_labels,
                                   k_sweep, loo_n, loo_mode, agg_fn,
                                   sample_group = NULL, progress_fn = NULL) {
  all_cts        <- colnames(niche_mat)
  unique_samples <- unique(sample_labels)
  n_s            <- length(unique_samples)
  n_iters        <- min(n_s, 20L)

  if (loo_mode == "group") {
    # sample_group: named vector, names = sample names, values = group label
    # (e.g. timepoint, from the mapped condition). Holds out loo_n samples
    # from *every* group each fold, instead of drawing loo_n/pct from the
    # pooled sample list, so a stratified structure (e.g. N timepoints x R
    # replicates) actually gets tested as intended.
    groups       <- sample_group[unique_samples]
    group_names  <- unique(groups)
    group_idx    <- lapply(group_names, function(g) unique_samples[groups == g])
    min_group_sz <- min(lengths(group_idx))
    n_hold_per_group <- max(1L, min(as.integer(loo_n), min_group_sz - 1L))
    if (n_hold_per_group < 1L || length(group_names) < 2L) {
      return(data.frame(k = k_sweep, stability = rep(0, length(k_sweep))))
    }
    loo_sets <- lapply(seq_len(n_iters), function(i) {
      set.seed(i)
      unlist(lapply(group_idx, function(idx) sample(idx, n_hold_per_group)))
    })
  } else {
    n_hold <- if (loo_mode == "count") min(as.integer(loo_n), n_s - 1L)
              else max(1L, floor(n_s * loo_n / 100))
    loo_sets <- lapply(seq_len(n_iters), function(i) {
      set.seed(i); sample(unique_samples, n_hold)
    })
  }

  scores <- numeric(length(k_sweep))
  for (ki in seq_along(k_sweep)) {
    if (!is.null(progress_fn)) progress_fn(ki, length(k_sweep), k_sweep[ki])
    k2 <- k_sweep[ki]
    iter_deltas <- numeric(n_iters)
    for (li in seq_along(loo_sets)) {
      held        <- loo_sets[[li]]
      train_mask  <- !sample_labels %in% held
      test_mask   <-  sample_labels %in% held
      train_data  <- niche_mat[train_mask, , drop = FALSE]
      if (nrow(train_data) < k2) { iter_deltas[li] <- NA_real_; next }
      km <- tryCatch(
        kmeans(train_data, centers = k2, nstart = 5L, iter.max = 100L),
        error = function(e) NULL)
      if (is.null(km)) { iter_deltas[li] <- NA_real_; next }
      tr_assign  <- km$cluster
      te_assign  <- if (sum(test_mask))
        .r_assign_to_centers(niche_mat[test_mask, , drop = FALSE], km$centers)
        else integer(0)
      all_assign           <- integer(nrow(niche_mat))
      all_assign[train_mask] <- tr_assign
      all_assign[test_mask]  <- te_assign
      ft <- compute_nh_freq(tr_assign, ctypes_vec[train_mask], k2, all_cts)
      ff <- compute_nh_freq(all_assign, ctypes_vec, k2, all_cts)
      iter_deltas[li] <- agg_fn(abs(ff - ft))
    }
    scores[ki] <- mean(iter_deltas, na.rm = TRUE)
  }
  gc()
  data.frame(k = k_sweep, stability = scores)
}

# ─── Unified niche-matrix builder ─────────────────────────────────────────────
# Returns a list with either a Python object ref or an R matrix for niche_mat,
# plus all metadata vectors.

build_niche_data <- function(spe_list, sample_names, ct_cols, k1) {
  # Global celltype universe
  all_cts <- sort(unique(unlist(lapply(sample_names, function(s) {
    ct_col <- ct_cols[[s]]
    sort(unique(as.character(colData(spe_list[[s]])[[ct_col]])))
  }))))

  ct_raw_list  <- vector("list", length(sample_names))
  ct_enc_list  <- vector("list", length(sample_names))  # 0-indexed global
  coords_list  <- vector("list", length(sample_names))
  n_cells_v    <- integer(length(sample_names))

  for (i in seq_along(sample_names)) {
    s       <- sample_names[i]
    spe     <- spe_list[[s]]
    ct_col  <- ct_cols[[s]]
    ctypes  <- as.character(colData(spe)[[ct_col]])
    ct_raw_list[[i]]  <- ctypes
    ct_enc_list[[i]]  <- match(ctypes, all_cts) - 1L   # 0-indexed
    coords_list[[i]]  <- spatialCoords(spe)
    n_cells_v[i]      <- ncol(spe)
  }

  cell_types_v  <- unlist(ct_raw_list)
  ct_encoded_v  <- unlist(ct_enc_list)
  sample_labels <- rep(sample_names, n_cells_v)
  samp_encoded  <- rep(seq_along(sample_names) - 1L, n_cells_v)  # 0-indexed

  if (PY_AVAILABLE) {
    # Stays as float32 numpy array in Python heap — no R copy until assignments
    py_nm <- reticulate::py$build_pooled_niche_matrix(
      coords_list, ct_enc_list, length(all_cts), as.integer(k1)
    )
    return(list(
      niche_mat     = py_nm,
      is_python     = TRUE,
      cell_types_v  = cell_types_v,
      ct_encoded_v  = ct_encoded_v,
      sample_labels = sample_labels,
      samp_encoded  = samp_encoded,
      all_cts       = all_cts,
      n_samples     = length(sample_names)
    ))
  }

  # ── R fallback ──
  mats <- lapply(seq_along(sample_names), function(i) {
    s  <- sample_names[i]
    m  <- .r_build_niche_matrix(spe_list[[s]], ct_cols[[s]], k1)
    missing_cts <- setdiff(all_cts, colnames(m))
    if (length(missing_cts)) {
      pad <- matrix(0, nrow(m), length(missing_cts),
                    dimnames = list(NULL, missing_cts))
      m   <- cbind(m, pad)
    }
    m[, all_cts, drop = FALSE]
  })
  niche_mat <- do.call(rbind, mats)
  colnames(niche_mat) <- all_cts
  rm(mats); gc()

  list(
    niche_mat     = niche_mat,
    is_python     = FALSE,
    cell_types_v  = cell_types_v,
    ct_encoded_v  = ct_encoded_v,
    sample_labels = sample_labels,
    samp_encoded  = samp_encoded,
    all_cts       = all_cts,
    n_samples     = length(sample_names)
  )
}

# ─── Other helpers ────────────────────────────────────────────────────────────
find_celltype_cols <- function(spe) {
  cols <- colnames(colData(spe))
  grep("annotation|celltype|cell_type|cluster|label", cols,
       ignore.case = TRUE, value = TRUE)
}

# ─── CSS ─────────────────────────────────────────────────────────────────────
APP_CSS <- "
body          { background:#f4f6fb; font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif; }
.panel-card   { background:#fff; border-radius:10px; padding:22px 26px;
                box-shadow:0 2px 8px rgba(0,0,0,.08); margin-bottom:18px; }
.panel-title  { font-size:1.05em; font-weight:700; color:#2d3a5c;
                border-left:4px solid #4e79a7; padding-left:10px; margin-bottom:14px; }
.step-badge   { display:inline-block; background:#4e79a7; color:#fff; border-radius:50%;
                width:24px; height:24px; text-align:center; line-height:24px;
                font-size:.82em; font-weight:700; margin-right:7px; }
.info-tag     { font-size:.78em; color:#6b7a99; margin-top:3px; }
hr.sep        { border-top:1px solid #e8ecf3; margin:16px 0; }
.shiny-notification { position:fixed; top:80%; left:40%; opacity:1;
                      height:auto; min-height:50px; width:360px; }
.nav-tabs>li>a        { color:#4e79a7; font-weight:600; }
.nav-tabs>li.active>a { color:#2d3a5c; border-top:3px solid #4e79a7; }
.btn-run      { background:#4e79a7; color:#fff; border:none; border-radius:6px;
                padding:8px 22px; font-weight:600; }
.btn-run:hover{ background:#3b618f; color:#fff; }
.btn-dl       { background:#59a14f; color:#fff; border:none; border-radius:6px;
                padding:8px 22px; font-weight:600; }
.btn-dl:hover { background:#468040; color:#fff; }
.optimal-k    { font-size:1.4em; font-weight:700; color:#e15759; }
.mem-bar      { font-size:.75em; color:#888; background:#f0f2f8;
                border-radius:4px; padding:3px 10px; display:inline-block; }
"

# ═══════════════════════════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════════════════════════
ui <- fluidPage(
  tags$head(
    tags$title("NeighborhoodR"),
    tags$style(HTML(APP_CSS))
  ),

  titlePanel(div(
    tags$span("NeighborhoodR",
              style = "font-weight:800; font-size:1.55em; color:#2d3a5c;"),
    tags$span(" — spatial neighborhood analysis",
              style = "color:#7a8ab0; font-size:.95em; margin-left:8px;")
  )),

  navbarPage(
    title = NULL, id = "main_tabs",

    # ── Tab 1: Upload ────────────────────────────────────────────────────────
    tabPanel("1 · Upload & Configure",
      fluidRow(
        column(4,
          div(class = "panel-card",
            div(class = "panel-title",
                tags$span(class = "step-badge", "1"), "Load Samples"),
            fileInput("rds_files", "Upload .rds SpatialExperiment files",
                      accept = ".rds", multiple = TRUE),
            div(class = "info-tag",
                "Multiple files accepted. Each must be a SpatialExperiment object."),
            hr(class = "sep"),
            div(class = "panel-title",
                tags$span(class = "step-badge", "2"), "KNN Parameters"),
            numericInput("k1", "K₁ (niche composition neighbours)", value = 10,
                         min = 3, max = 500),
            div(class = "info-tag",
                "Spatial nearest-neighbours used to build each cell's niche composition vector.")
          ),
          div(class = "panel-card",
            div(class = "panel-title", "Memory"),
            textOutput("mem_display"),
            br(),
            div(class = "info-tag",
                if (PY_AVAILABLE)
                  "Niche matrix stored as float32 in Python heap (~50 % less RAM than R default)."
                else
                  "R backend: matrix stored as float64. Install sklearn/scipy for memory savings."),
            hr(class = "sep"),
            actionButton("reset_session", "⟳ New analysis (clear memory)",
                         style = "background:none; border:none; color:#aaa; font-size:.8em;
                                  padding:2px 0; cursor:pointer; text-decoration:underline;",
                         title = "Null all loaded objects and reload the session")
          )
        ),
        column(8,
          div(class = "panel-card",
            div(class = "panel-title",
                tags$span(class = "step-badge", "3"),
                "Cell-type Column — per Sample"),
            div(class = "info-tag",
                "Select the metadata column containing cell-type/cluster annotations. ★ marks auto-detected columns."),
            br(),
            uiOutput("celltype_selectors"),
            conditionalPanel("output.samples_loaded",
              hr(class = "sep"),
              div(class = "panel-title",
                  tags$span(class = "step-badge", "4"),
                  "Condition Column (optional)"),
              uiOutput("condition_selector"),
              div(class = "info-tag",
                  "Select a shared metadata column, or leave as (none) and assign labels below."),
              hr(class = "sep"),
              div(class = "panel-title",
                  tags$span(class = "step-badge", "5"),
                  "Sample → Condition Map"),
              div(class = "info-tag",
                  "Assign a condition label to each sample. Used only when no metadata column is selected above; ignored if both are set."),
              br(),
              uiOutput("condition_map_ui")
            )
          )
        )
      )
    ),

    # ── Tab 2: Stability sweep ───────────────────────────────────────────────
    tabPanel("2 · Neighborhood Optimisation",
      fluidRow(
        column(4,
          div(class = "panel-card",
            div(class = "panel-title",
                tags$span(class = "step-badge", "A"), "K₂ Sweep Range"),
            numericInput("k2_min", "K₂ minimum", value = 3, min = 2),
            uiOutput("k2_max_ui"),
            div(class = "info-tag",
                "Sweep from K₂ min to K₂ max — defaults to [3, N_celltypes − 1]."),
            hr(class = "sep"),
            div(class = "panel-title",
                tags$span(class = "step-badge", "B"), "Leave-one-out Parameters"),
            radioButtons("loo_mode", "Hold-out unit",
                         choices = c("Number of samples" = "count",
                                     "Percentage of samples" = "pct",
                                     "Per mapped condition" = "group"),
                         inline = TRUE),
            conditionalPanel("input.loo_mode == 'count'",
              numericInput("loo_n_count", "Samples to hold out", value = 1, min = 1)
            ),
            conditionalPanel("input.loo_mode == 'pct'",
              numericInput("loo_n_pct", "% of samples to hold out",
                           value = 20, min = 1, max = 90)
            ),
            conditionalPanel("input.loo_mode == 'group'",
              numericInput("loo_n_group", "Samples to hold out per group", value = 1, min = 1),
              div(class = "info-tag",
                  "Groups are the Condition Column / Sample → Condition Map set in Tab 1 ",
                  "(e.g. timepoint) — every fold holds out this many samples from ",
                  "every group, instead of drawing from the pooled sample list.")
            ),
            selectInput("agg_fn", "Stability metric",
                        choices = c("Median absolute delta" = "median",
                                    "Mean absolute delta"   = "mean")),
            hr(class = "sep"),
            actionButton("run_sweep", "Run Stability Sweep",
                         class = "btn-run", width = "100%"),
            br(), br(),
            div(class = "info-tag",
                "Optimal K₂ minimises the curve — neighborhoods are most stable when adding held-out samples changes their celltype composition the least.")
          ),
          div(class = "panel-card",
            div(class = "panel-title", "Memory After Sweep"),
            textOutput("mem_after_sweep")
          )
        ),
        column(8,
          div(class = "panel-card",
            div(class = "panel-title", "Stability Curve"),
            plotOutput("stability_plot", height = "380px") %>%
              withSpinner(color = "#4e79a7"),
            br(),
            div(style = "text-align:center;", uiOutput("optimal_k_display"))
          ),
          div(class = "panel-card",
            div(class = "panel-title", "Sweep Table"),
            DTOutput("stability_table") %>% withSpinner(color = "#4e79a7")
          )
        )
      )
    ),

    # ── Tab 3: Assign & Download ─────────────────────────────────────────────
    tabPanel("3 · Assign & Download",
      fluidRow(
        column(4,
          div(class = "panel-card",
            div(class = "panel-title",
                tags$span(class = "step-badge", "!"), "Final K₂"),
            uiOutput("final_k2_ui"),
            div(class = "info-tag",
                "Pre-filled from the sweep optimum. Override if desired."),
            hr(class = "sep"),
            actionButton("run_assign", "Assign Neighborhoods",
                         class = "btn-run", width = "100%"),
            br(), br(),
            div(class = "info-tag",
                "Already assigned? Update condition labels in Tab 1, then click below to refresh plots without re-running."),
            actionButton("update_condition", "Update Condition Map",
                         style = "background:#5a9f6e; color:#fff; border:none; border-radius:6px;
                                  padding:6px 12px; font-weight:600; font-size:.85em;
                                  cursor:pointer; width:100%; margin-top:4px;"),
            br(), br(),
            uiOutput("download_ui")
          )
        ),
        column(8,
          div(class = "panel-card",
            div(class = "panel-title", "Assignment Summary"),
            DTOutput("assignment_summary") %>% withSpinner(color = "#4e79a7")
          )
        )
      )
    ),

    # ── Tab 4: Visualisations ────────────────────────────────────────────────
    tabPanel("4 · Visualisations",
      fluidRow(
        column(3,
          div(class = "panel-card",
            div(class = "panel-title", "Plot Controls"),
            selectInput("viz_sample",  "Sample", choices = NULL),
            selectInput("viz_colour", "Colour cells by",
                        choices = c("Neighborhood" = "neighborhood",
                                    "Cell type"     = "celltype")),
            numericInput("pt_size", "Point size",  value = 1.5,
                         min = 0.2, max = 8,  step = 0.1),
            sliderInput("pt_alpha",  "Opacity",
                        min = 0.1, max = 1,   value = 0.8, step = 0.05)
          )
        ),
        column(9,
          div(class = "panel-card",
            div(class = "panel-title", "Spatial Projection"),
            plotOutput("spatial_plot", height = "500px") %>%
              withSpinner(color = "#4e79a7")
          )
        )
      ),
      fluidRow(
        column(12,
          div(class = "panel-card",
            div(class = "panel-title",
                "Celltype Composition per Neighborhood (all samples)"),
            plotOutput("barplot_nh", height = "420px") %>%
              withSpinner(color = "#4e79a7")
          )
        )
      ),
      conditionalPanel("output.has_condition",
        fluidRow(
          column(12,
            div(class = "panel-card",
              div(class = "panel-title", "Neighborhood Abundance by Condition"),
              plotOutput("condition_plot", height = "420px") %>%
                withSpinner(color = "#4e79a7"),
              br(),
              div(class = "panel-title", "Kruskal-Wallis Statistics"),
              DTOutput("condition_stats") %>% withSpinner(color = "#4e79a7")
            )
          )
        )
      )
    ),

    # ── Tab 5: Reproducibility ───────────────────────────────────────────────
    tabPanel("5 · Reproducibility",
      fluidRow(
        column(6,
          div(class = "panel-card",
            div(class = "panel-title", "Session & Parameters"),
            verbatimTextOutput("session_info_out")
          )
        ),
        column(6,
          div(class = "panel-card",
            div(class = "panel-title", "Replay Script"),
            div(class = "info-tag",
                "Stand-alone R script reproducing this analysis from the same input files."),
            br(),
            downloadButton("download_script",   "Download Replay Script", class = "btn-dl"),
            br(), br(),
            downloadButton("download_prov_json","Download Provenance JSON", class = "btn-dl")
          ),
          div(class = "panel-card",
            div(class = "panel-title", "Provenance JSON Preview"),
            verbatimTextOutput("prov_preview")
          )
        )
      )
    )
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# Server
# ═══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  prov <- if (exists("ProvenanceTracker")) {
    ProvenanceTracker$new("NeighborhoodR", session = session, output_dir = tempdir())
  } else NULL

  rv <- reactiveValues(
    spe_list      = NULL,
    sample_names  = NULL,
    # niche_mat is either a Python numpy object (PY_AVAILABLE) or an R matrix
    niche_mat     = NULL,
    is_python     = FALSE,
    cell_types_v  = NULL,
    ct_encoded_v  = NULL,
    sample_labels = NULL,
    samp_encoded  = NULL,
    all_cts       = NULL,
    n_samples     = NULL,
    sweep_res     = NULL,
    optimal_k2    = NULL,
    assignments   = NULL,
    joint_spe     = NULL,
    spe_list_full = NULL,  # full per-sample SPEs (markers intact) for export
    has_condition = FALSE,
    condition_col = NULL,
    condition_map = NULL,   # named list: sample -> label (from manual map UI)
    mem_sweep     = NULL
  )

  # ── Memory display ─────────────────────────────────────────────────────────
  mem_tick <- reactiveTimer(8000)  # refresh every 8 s

  output$mem_display <- renderText({
    mem_tick()
    mem_label()
  })
  output$mem_after_sweep <- renderText({
    rv$mem_sweep %||% "—"
  })

  # ── 1. Load RDS files ──────────────────────────────────────────────────────
  observeEvent(input$rds_files, {
    req(input$rds_files)
    fi <- input$rds_files

    # Clear any previously computed niche matrix from Python heap
    if (rv$is_python) {
      rv$niche_mat <- NULL
      if (PY_AVAILABLE) reticulate::py_run_string("_gc.collect()")
    }
    gc()

    withProgress(message = "Loading .rds files…", value = 0, {
      spe_list <- list()
      for (i in seq_len(nrow(fi))) {
        incProgress(1 / nrow(fi), detail = fi$name[i])
        obj <- tryCatch(readRDS(fi$datapath[i]), error = function(e) NULL)
        if (is.null(obj) || !is(obj, "SpatialExperiment")) {
          showNotification(
            paste(fi$name[i], "— not a SpatialExperiment, skipped."),
            type = "warning")
          next
        }
        # Keep the full object (markers intact) here — it's what gets
        # exported for downstream analysis/visualisation outside this app.
        # Marker-free copies are built on demand at the concatenation step
        # (see run_assign), since per-sample marker panels can differ and
        # cbind() requires identical row counts across objects.
        sname <- tools::file_path_sans_ext(fi$name[i])
        spe_list[[sname]] <- obj
        if (!is.null(prov)) prov$register_input(fi[i, ], input_id = "rds_files")
      }
    })

    if (!length(spe_list)) {
      showNotification("No valid SpatialExperiment objects loaded.", type = "error")
      return()
    }
    rv$spe_list     <- spe_list
    rv$sample_names <- names(spe_list)
    updateSelectInput(session, "viz_sample", choices = rv$sample_names)
  })

  # ── 2. Dynamic per-sample celltype selectors ───────────────────────────────
  output$samples_loaded <- reactive({ !is.null(rv$spe_list) })
  outputOptions(output, "samples_loaded", suspendWhenHidden = FALSE)

  output$celltype_selectors <- renderUI({
    req(rv$spe_list)
    lapply(rv$sample_names, function(sname) {
      all_cols  <- colnames(colData(rv$spe_list[[sname]]))
      auto_cols <- find_celltype_cols(rv$spe_list[[sname]])
      div(style = "margin-bottom:12px;",
        selectInput(
          inputId  = paste0("ct_col_", make.names(sname)),
          label    = paste0(sname, " — cell-type column"),
          choices  = setNames(all_cols,
                              ifelse(all_cols %in% auto_cols,
                                     paste0("★ ", all_cols), all_cols)),
          selected = auto_cols[1] %||% all_cols[1]
        )
      )
    })
  })

  output$condition_selector <- renderUI({
    req(rv$spe_list)
    shared <- Reduce(intersect,
                     lapply(rv$spe_list, function(s) colnames(colData(s))))
    tagList(
      selectInput("condition_col", "Condition column",
                  choices  = c("(none)" = "", shared),
                  selected = grep("condition", shared, ignore.case = TRUE,
                                  value = TRUE)[1] %||% "")
    )
  })

  output$condition_map_ui <- renderUI({
    req(rv$spe_list)
    lapply(rv$sample_names, function(sname) {
      existing <- rv$condition_map[[sname]] %||% ""
      div(style = "margin-bottom:8px;",
        textInput(
          inputId     = paste0("cond_map_", make.names(sname)),
          label       = sname,
          value       = existing,
          placeholder = "e.g. treated, control, …"
        )
      )
    })
  })

  # ── 3. K₂ max default based on celltype count ──────────────────────────────
  output$k2_max_ui <- renderUI({
    n_def <- if (!is.null(rv$all_cts)) length(rv$all_cts) - 1L else 10L
    numericInput("k2_max", "K₂ maximum",
                 value = max(n_def, (input$k2_min %||% 3L) + 1L),
                 min   = (input$k2_min %||% 2L) + 1L)
  })

  # ── 4. Helper: gather ct_cols named vector ─────────────────────────────────
  get_ct_cols <- reactive({
    req(rv$spe_list)
    vapply(rv$sample_names, function(s) {
      id  <- paste0("ct_col_", make.names(s))
      input[[id]] %||% find_celltype_cols(rv$spe_list[[s]])[1]
    }, character(1))
  })

  # ── 5. Stability sweep ────────────────────────────────────────────────────
  observeEvent(input$run_sweep, {
    req(rv$spe_list, input$k1, input$k2_min, input$k2_max)

    ct_cols <- isolate(get_ct_cols())

    # ── Build pooled niche matrix ──
    withProgress(message = "Building niche matrix…", value = 0.1, {
      nm_data <- tryCatch(
        build_niche_data(rv$spe_list, rv$sample_names, ct_cols, input$k1),
        error = function(e) {
          showNotification(paste("Niche matrix error:", e$message), type = "error")
          NULL
        }
      )
    })
    req(nm_data)

    rv$niche_mat     <- nm_data$niche_mat
    rv$is_python     <- nm_data$is_python
    rv$cell_types_v  <- nm_data$cell_types_v
    rv$ct_encoded_v  <- nm_data$ct_encoded_v
    rv$sample_labels <- nm_data$sample_labels
    rv$samp_encoded  <- nm_data$samp_encoded
    rv$all_cts       <- nm_data$all_cts
    rv$n_samples     <- nm_data$n_samples
    rm(nm_data); gc()

    # Update K₂ max if needed
    updateNumericInput(session, "k2_max",
                       value = max(input$k2_max, length(rv$all_cts) - 1L))

    k_sweep  <- seq(as.integer(input$k2_min), as.integer(input$k2_max))
    loo_n    <- switch(input$loo_mode,
                        count = input$loo_n_count,
                        pct   = input$loo_n_pct,
                        group = input$loo_n_group)

    # ── Group mode: resolve one label per sample from the Condition Column /
    # Sample -> Condition Map set in Tab 1, same source as .cell_condition().
    sample_group_v <- NULL
    if (input$loo_mode == "group") {
      cond_col <- input$condition_col
      groups <- if (!is.null(cond_col) && nchar(cond_col)) {
        vapply(rv$sample_names, function(s) {
          v <- as.character(colData(rv$spe_list[[s]])[[cond_col]])
          v <- v[!is.na(v) & nchar(v)]
          if (length(v)) v[1] else NA_character_
        }, character(1))
      } else {
        vapply(rv$sample_names, function(s) {
          v <- trimws(input[[paste0("cond_map_", make.names(s))]] %||% "")
          if (nchar(v)) v else NA_character_
        }, character(1))
      }
      names(groups) <- rv$sample_names
      if (anyNA(groups)) {
        showNotification(
          "Per mapped condition needs every sample assigned — set the Condition Column or fill in every Sample -> Condition Map row in Tab 1.",
          type = "error")
        return()
      }
      if (length(unique(groups)) < 2L) {
        showNotification(
          "Per mapped condition needs at least 2 distinct condition values across your samples.",
          type = "error")
        return()
      }
      sample_group_v <- groups
    }

    # ── Run sweep ──
    withProgress(message = "Running LOO stability sweep…", value = 0, {
      if (rv$is_python) {
        sample_group_enc <- if (!is.null(sample_group_v))
          as.integer(match(sample_group_v[rv$sample_names], unique(sample_group_v)) - 1L)
        else NULL
        res_py <- reticulate::py$loo_stability_sweep(
          rv$niche_mat,
          as.integer(rv$ct_encoded_v),
          as.integer(rv$samp_encoded),
          as.integer(rv$n_samples),
          as.integer(k_sweep),
          loo_n, input$loo_mode, input$agg_fn,
          as.integer(length(rv$all_cts)),
          sample_group_r = sample_group_enc
        )
        sweep_res <- data.frame(
          k         = as.integer(res_py$k),
          stability = as.numeric(res_py$stability)
        )
        reticulate::py_run_string("_gc.collect()")
      } else {
        agg_fn <- if (input$agg_fn == "median") median else mean
        sweep_res <- .r_loo_stability_sweep(
          rv$niche_mat, rv$cell_types_v, rv$sample_labels,
          k_sweep, loo_n, input$loo_mode, agg_fn,
          sample_group = sample_group_v,
          progress_fn = function(ki, total, k_val)
            incProgress(1 / total, detail = paste0("k = ", k_val))
        )
      }
    })

    rv$sweep_res  <- sweep_res
    best_idx      <- which.min(sweep_res$stability)
    rv$optimal_k2 <- if (length(best_idx)) sweep_res$k[best_idx] else sweep_res$k[1L]
    rv$mem_sweep  <- mem_label()
    updateNumericInput(session, "final_k2_override", value = rv$optimal_k2)

    if (!is.null(prov)) {
      prov$capture_parameters(reactiveValuesToList(input))
      prov$analysis_started()
    }
    showNotification(paste0("Sweep complete. Optimal K₂ = ", rv$optimal_k2),
                     type = "message")
  })

  # ── 6. Stability outputs ───────────────────────────────────────────────────
  output$stability_plot <- renderPlot({
    req(rv$sweep_res)
    df  <- rv$sweep_res
    opt <- rv$optimal_k2
    ggplot(df, aes(k, stability)) +
      geom_line(colour = "#4e79a7", linewidth = 1.1) +
      geom_point(colour = "#4e79a7", size = 2.5) +
      geom_vline(xintercept = opt, colour = "#e15759",
                 linetype = "dashed", linewidth = 1) +
      annotate("text", x = opt, y = max(df$stability, na.rm = TRUE),
               label = paste0(" K₂ = ", opt), colour = "#e15759",
               hjust = -0.1, fontface = "bold") +
      labs(x = "K₂ (number of neighborhoods)",
           y = paste(ifelse(input$agg_fn == "median", "Median", "Mean"),
                     "absolute Δ (neighborhood × celltype frequency)"),
           title = "LOO Stability Curve") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$stability_table <- renderDT({
    req(rv$sweep_res)
    df <- rv$sweep_res
    df$stability <- round(df$stability, 6)
    opt          <- rv$optimal_k2
    df$optimal   <- if (length(opt)) df$k == opt else FALSE
    datatable(df, rownames = FALSE,
              options = list(pageLength = 10, dom = "tp")) %>%
      formatStyle("optimal", target = "row",
                  backgroundColor = styleEqual(TRUE, "#fff3cd"))
  })

  output$optimal_k_display <- renderUI({
    req(rv$optimal_k2)
    div(
      tags$span("Optimal K₂: ", style = "color:#444; font-size:1.1em;"),
      tags$span(class = "optimal-k", rv$optimal_k2),
      tags$span(" neighborhoods", style = "color:#444; font-size:1.1em;")
    )
  })

  # ── 7. Assign neighborhoods ──────────────────────────────────────────────
  output$final_k2_ui <- renderUI({
    numericInput("final_k2_override", "K₂ for final assignment",
                 value = rv$optimal_k2 %||% 6L, min = 2L)
  })

  observeEvent(input$run_assign, {
    req(rv$niche_mat, input$final_k2_override)
    k2 <- as.integer(input$final_k2_override)

    withProgress(message = "Assigning neighborhoods…", value = 0.15, {
      if (rv$is_python) {
        km_res      <- reticulate::py$final_kmeans(rv$niche_mat, k2)
        assignments <- as.integer(unlist(km_res$labels))
        reticulate::py_run_string("_gc.collect()")
      } else {
        set.seed(42L)
        km          <- kmeans(rv$niche_mat, centers = k2, nstart = 25L, iter.max = 300L)
        assignments <- km$cluster
        rm(km); gc()
      }

      incProgress(0.5, detail = "Writing to SPE objects")
      spe_out <- rv$spe_list
      ct_cols <- isolate(get_ct_cols())
      ptr     <- 1L
      for (sname in rv$sample_names) {
        n_cells <- ncol(spe_out[[sname]])
        colData(spe_out[[sname]])$neighborhood <-
          paste0("N", assignments[ptr:(ptr + n_cells - 1L)])
        colData(spe_out[[sname]])$sample <- rep(sname, n_cells)
        # Per-sample cell-type columns aren't guaranteed to share a name
        # across samples (e.g. cluster-label columns embed the per-sample
        # cluster count in their name), so normalise into a common field
        # for cross-sample use.
        colData(spe_out[[sname]])$celltype <-
          as.character(colData(spe_out[[sname]])[[ct_cols[[sname]]]])
        ptr <- ptr + n_cells
      }

      incProgress(0.25, detail = "Finalising")
      rv$assignments <- assignments

      # ── Condition: metadata column takes precedence over manual map ──
      cond_col <- input$condition_col
      cmap <- setNames(
        vapply(rv$sample_names, function(s) {
          v <- input[[paste0("cond_map_", make.names(s))]] %||% ""
          trimws(v)
        }, character(1)),
        rv$sample_names
      )
      cmap_any <- any(nchar(cmap) > 0)

      if (!is.null(cond_col) && nchar(cond_col)) {
        rv$condition_col <- cond_col
        rv$condition_map <- NULL
        rv$has_condition <- TRUE
      } else if (cmap_any) {
        rv$condition_col <- NULL
        rv$condition_map <- as.list(cmap)
        for (sname in rv$sample_names) {
          lbl <- cmap[[sname]]
          colData(spe_out[[sname]])$condition_map <-
            rep(if (nchar(lbl)) lbl else NA_character_, ncol(spe_out[[sname]]))
        }
        rv$has_condition <- TRUE
      }

      # spe_out keeps full marker data — retained as-is for downstream
      # export (see download_bundle). It's not used for anything internal
      # to this app, so it never needs to be concatenated as-is.
      rv$spe_list_full <- spe_out

      # For internal use (plots, summary table, joint object) build
      # marker-free copies just for this concatenation: cbind() requires
      # identical row counts AND identical colData column names across
      # objects. Marker panels can differ per sample (zero-variance markers
      # dropped independently upstream), and original colData schemas can
      # differ too (e.g. per-sample cluster-label columns with different
      # names) — none of which this app's own logic reads. Keep only the
      # uniform fields every sample now has.
      light_cols <- c("sample", "neighborhood", "celltype",
                       if (!is.null(cond_col) && nchar(cond_col)) cond_col,
                       "condition_map")
      spe_light <- lapply(spe_out, function(s) {
        cd   <- colData(s)
        keep <- intersect(light_cols, colnames(cd))
        SpatialExperiment::SpatialExperiment(
          assays        = list(placeholder = matrix(numeric(0), nrow = 0, ncol = ncol(s))),
          colData       = cd[, keep, drop = FALSE],
          spatialCoords = spatialCoords(s)
        )
      })
      rv$joint_spe <- do.call(cbind, unname(spe_light))
      rm(spe_out, spe_light, assignments); gc()
    })

    if (!is.null(prov)) {
      prov$capture_parameters(reactiveValuesToList(input))
      prov$analysis_completed()
    }
    showNotification("Neighborhood assignment complete.", type = "message")
    updateTabsetPanel(session, "main_tabs", selected = "4 · Visualisations")
  })

  output$has_condition <- reactive({ isTRUE(rv$has_condition) })
  outputOptions(output, "has_condition", suspendWhenHidden = FALSE)

  # ── Shared summary reactives ────────────────────────────────────────────
  # Computed once here and reused by both the on-screen outputs below and
  # the zip downloads, so a download always matches what the GUI shows.
  assignment_summary_df <- reactive({
    if (is.null(rv$joint_spe)) return(NULL)
    cd   <- as.data.frame(colData(rv$joint_spe)[, c("sample", "neighborhood"), drop = FALSE])
    rows <- lapply(rv$sample_names, function(sname) {
      tab <- table(cd$neighborhood[cd$sample == sname])
      data.frame(Sample = sname, Neighborhood = names(tab), N_cells = as.integer(tab))
    })
    do.call(rbind, rows)
  })

  barplot_nh_plot <- reactive({
    if (is.null(rv$joint_spe) || is.null(rv$all_cts)) return(NULL)
    rows <- lapply(rv$sample_names, function(sname) {
      spe <- rv$joint_spe[, colData(rv$joint_spe)$sample == sname]
      if (!"celltype" %in% colnames(colData(spe))) return(NULL)
      data.frame(neighborhood = colData(spe)$neighborhood,
                 celltype      = as.character(colData(spe)$celltype),
                 sample        = sname)
    })
    df <- do.call(rbind, Filter(Negate(is.null), rows))
    if (is.null(df) || !nrow(df)) return(NULL)

    df_freq <- df %>%
      count(neighborhood, celltype) %>%
      group_by(neighborhood) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup()

    all_cts_plot <- sort(unique(df_freq$celltype))
    pal_vec      <- setNames(pal(length(all_cts_plot)), all_cts_plot)

    ggplot(df_freq, aes(neighborhood, prop, fill = celltype)) +
      geom_col(position = "stack", width = 0.75) +
      scale_fill_manual(values = pal_vec) +
      scale_y_continuous(labels = percent_format()) +
      labs(x = "Neighborhood", y = "Proportion", fill = "Cell type",
           title = "Celltype composition per neighborhood") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            plot.title  = element_text(face = "bold"))
  })

  condition_stats_df <- reactive({
    if (is.null(rv$joint_spe) || !isTRUE(rv$has_condition)) return(NULL)
    rows <- lapply(rv$sample_names, function(sname) {
      spe  <- rv$joint_spe[, colData(rv$joint_spe)$sample == sname]
      cond <- .cell_condition(sname, spe)
      if (is.null(cond)) return(NULL)
      data.frame(neighborhood = colData(spe)$neighborhood,
                 condition     = cond,
                 sample        = sname)
    })
    df <- do.call(rbind, Filter(Negate(is.null), rows))
    if (is.null(df) || !nrow(df)) return(NULL)

    df_freq <- df %>%
      count(sample, condition, neighborhood) %>%
      group_by(sample) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup()

    df_freq %>%
      group_by(neighborhood) %>%
      summarise(
        kw_statistic = tryCatch(
          round(kruskal.test(prop ~ condition, data = cur_data())$statistic, 4),
          error = function(e) NA_real_),
        p_value = tryCatch(
          signif(kruskal.test(prop ~ condition, data = cur_data())$p.value, 4),
          error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      mutate(p_adj = signif(p.adjust(p_value, method = "BH"), 4),
             significant = p_adj < 0.05)
  })

  output$assignment_summary <- renderDT({
    df <- assignment_summary_df()
    req(df)
    datatable(df, rownames = FALSE, options = list(pageLength = 15, dom = "ftp"))
  })

  # estimated in-memory size (bytes) of the full-marker export, since that's
  # what actually gets written in the under-threshold branch below
  .joint_size <- reactive({
    req(rv$spe_list_full)
    tryCatch(as.numeric(object.size(rv$spe_list_full)), error = function(e) 0)
  })

  output$download_ui <- renderUI({
    req(rv$joint_spe, rv$spe_list_full)
    sz   <- .joint_size()
    big  <- sz >= 5e9
    lbl  <- if (big)
      paste0("Download bundle (.zip, ~", round(sz / 1e9, 1), " GB coldata + coordinates + summaries)")
    else
      paste0("Download bundle (.zip, ~", round(sz / 1e6), " MB SPE objects + coldata + coordinates + summaries)")
    tagList(
      br(),
      div(class = "info-tag",
          if (big)
            paste0("Object exceeds 5 GB — bundle skips the per-sample SPE .rds objects (markers), but ",
                   "includes per-sample coldata tables (x, y, neighborhood, celltype, condition), ",
                   "per-sample/per-neighborhood coordinate files for log-odds and PCF analysis, the ",
                   "assignment summary, celltype composition plot, and Kruskal-Wallis condition statistics.")
          else
            paste0("Under 5 GB — bundle contains one full SpatialExperiment (.rds) per sample (markers ",
                   "intact), per-sample coldata tables (x, y, neighborhood, celltype, condition), ",
                   "per-sample/per-neighborhood coordinate files for log-odds and PCF analysis, the ",
                   "assignment summary, celltype composition plot, and Kruskal-Wallis condition statistics.")),
      br(),
      downloadButton("download_bundle", lbl, class = "btn-dl")
    )
  })

  # ── Download bundle ─────────────────────────────────────────────────────
  # One .zip covering everything: per-sample SPE objects (markers intact,
  # skipped only when the full-marker export would exceed 5 GB — a server-
  # side memory/disk concern while building the zip, not a browser download
  # limit), per-sample coldata tables, per-sample/per-neighborhood coordinate
  # files for log-odds/PCF tools, and the summary table/stats/plot. Kept as a
  # single download rather than split by content — uploads are what's capped
  # (shiny.maxRequestSize), not downloads, so there's no reason to make the
  # user grab two separate zips for one analysis.
  output$download_bundle <- downloadHandler(
    filename = function() paste0("neighborhoodR_", Sys.Date(), ".zip"),
    content  = function(file) {
      req(rv$joint_spe)
      tmp <- tempfile("neighborhoodR_export_")
      dir.create(tmp)
      on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

      include_spe <- .joint_size() < 5e9

      for (sname in rv$sample_names) {
        spe  <- rv$joint_spe[, colData(rv$joint_spe)$sample == sname]
        xy   <- as.data.frame(spatialCoords(spe))
        colnames(xy) <- c("x", "y")
        ct   <- if ("celltype" %in% colnames(colData(spe)))
          as.character(colData(spe)$celltype) else NA_character_
        nbhd <- as.character(colData(spe)$neighborhood)
        cond <- .cell_condition(sname, spe)
        if (is.null(cond)) cond <- NA_character_

        # ── Full per-sample SpatialExperiment (markers intact) ──
        # One .rds per sample rather than cbind()'d — per-sample marker
        # panels aren't guaranteed to match, and cbind() requires identical
        # row counts and colData column names across objects.
        if (include_spe)
          saveRDS(rv$spe_list_full[[sname]], file.path(tmp, paste0(sname, ".rds")))

        # ── Per-sample coldata (x, y, neighborhood, celltype, condition) ──
        # lets this sample's neighborhood assignment be mapped straight back
        # onto the original object without opening the .rds in R.
        write.csv(
          data.frame(x = xy$x, y = xy$y, neighborhood = nbhd,
                     celltype = ct, condition = cond, stringsAsFactors = FALSE),
          file.path(tmp, paste0(sname, "_coldata.csv")), row.names = FALSE)

        # ── Per-sample, per-neighborhood coordinate files ──
        # One CSV per (sample, neighborhood), in two formats:
        #  - logodds/: x, y, cluster        (pairwise_logOdds() input)
        #  - pcf/:     Sample Name, Cell X/Y Position, Tissue Category, Phenotype
        #              (extract_data()/pcf() input, Vectra cell_seg_data format)
        # Restricting each file to one neighborhood lets these tools be run
        # on a single niche within a sample, rather than the whole tissue.
        lo_dir  <- file.path(tmp, sname, "logodds")
        pcf_dir <- file.path(tmp, sname, "pcf")
        dir.create(lo_dir,  recursive = TRUE)
        dir.create(pcf_dir, recursive = TRUE)

        for (n in sort(unique(nbhd))) {
          mask <- nbhd == n
          tag  <- paste0(sname, "_", n)

          write.csv(
            data.frame(x = xy$x[mask], y = xy$y[mask], cluster = ct[mask],
                       stringsAsFactors = FALSE),
            file.path(lo_dir, paste0(tag, ".csv")), row.names = FALSE)

          write.csv(
            data.frame(
              `Sample Name`     = tag,
              `Cell X Position` = xy$x[mask],
              `Cell Y Position` = xy$y[mask],
              `Tissue Category` = "All",
              Phenotype         = ct[mask],
              check.names = FALSE, stringsAsFactors = FALSE
            ),
            file.path(pcf_dir, paste0(tag, ".csv")), row.names = FALSE)
        }
      }

      # ── Assignment summary, condition stats, and composition plot ──
      summ <- assignment_summary_df()
      if (!is.null(summ))
        write.csv(summ, file.path(tmp, "assignment_summary.csv"), row.names = FALSE)

      stats <- condition_stats_df()
      if (!is.null(stats))
        write.csv(stats, file.path(tmp, "condition_stats_kruskal_wallis.csv"), row.names = FALSE)

      plt <- barplot_nh_plot()
      if (!is.null(plt))
        ggsave(file.path(tmp, "celltype_composition_per_neighborhood.png"),
               plot = plt, width = 10, height = 6, dpi = 150)

      zip::zip(zipfile = file, files = dir(tmp), root = tmp, compression_level = 1)
    }
  )

  # ── 8. Visualisations ─────────────────────────────────────────────────────

  output$spatial_plot <- renderPlot({
    req(rv$joint_spe, input$viz_sample)
    sname <- input$viz_sample
    spe   <- rv$joint_spe[, colData(rv$joint_spe)$sample == sname]
    df    <- as.data.frame(spatialCoords(spe))
    colnames(df) <- c("x", "y")

    if (input$viz_colour == "neighborhood") {
      df$colour <- colData(spe)$neighborhood
    } else {
      df$colour <- if ("celltype" %in% colnames(colData(spe)))
        as.character(colData(spe)$celltype) else "unknown"
    }
    lvls    <- sort(unique(df$colour))
    pal_vec <- setNames(pal(length(lvls)), lvls)

    ggplot(df, aes(x, y, colour = colour)) +
      geom_point(size = input$pt_size, alpha = input$pt_alpha) +
      scale_colour_manual(values = pal_vec) +
      labs(title = paste0(sname, " — ", input$viz_colour),
           colour = input$viz_colour, x = "x", y = "y") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right", plot.title = element_text(face = "bold"))
  })

  output$barplot_nh <- renderPlot({
    p <- barplot_nh_plot()
    req(p)
    p
  })

  # Helper: extract per-cell condition vector for a sample
  .cell_condition <- function(sname, spe) {
    if (!is.null(rv$condition_col) && nchar(rv$condition_col) &&
        rv$condition_col %in% colnames(colData(spe)))
      return(as.character(colData(spe)[[rv$condition_col]]))
    if ("condition_map" %in% colnames(colData(spe)))
      return(as.character(colData(spe)$condition_map))
    NULL
  }

  output$condition_plot <- renderPlot({
    req(rv$joint_spe, rv$has_condition)
    cond_label <- rv$condition_col %||% "condition (mapped)"
    rows <- lapply(rv$sample_names, function(sname) {
      spe  <- rv$joint_spe[, colData(rv$joint_spe)$sample == sname]
      cond <- .cell_condition(sname, spe)
      if (is.null(cond)) return(NULL)
      data.frame(neighborhood = colData(spe)$neighborhood,
                 condition     = cond,
                 sample        = sname)
    })
    df <- do.call(rbind, Filter(Negate(is.null), rows))
    req(df)

    df_freq <- df %>%
      count(sample, condition, neighborhood) %>%
      group_by(sample) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup()

    conditions <- sort(unique(df_freq$condition))
    pal_vec    <- setNames(pal(length(conditions)), conditions)
    # Fix condition to a shared factor level set so position_dodge() lines up
    # identically between the boxplot layer and the singleton-point layer
    # below (each layer would otherwise compute dodge offsets from whatever
    # subset of conditions it happens to contain at each neighborhood).
    df_freq$condition <- factor(df_freq$condition, levels = conditions)

    # A boxplot needs >1 value to show anything but a flat, uncoloured line
    # (min = median = max collapses the box to zero height, so its fill
    # colour — the only thing the legend encodes — never renders). For any
    # (neighborhood, condition) backed by exactly one sample, overlay a
    # large coloured point at that value instead, so it's still readable
    # from the legend. Groups with >1 sample are untouched.
    df_freq <- df_freq %>%
      group_by(neighborhood, condition) %>%
      mutate(n_samples = n()) %>%
      ungroup()
    singles <- df_freq %>% filter(n_samples == 1)

    p <- ggplot(df_freq, aes(neighborhood, prop, fill = condition)) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.size = 1.2) +
      scale_fill_manual(values = pal_vec) +
      scale_y_continuous(labels = percent_format()) +
      labs(x = "Neighborhood", y = "Proportion of sample cells",
           fill = "Condition",
           title = paste0("Neighborhood abundance by ", cond_label)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            plot.title  = element_text(face = "bold"))

    if (nrow(singles) > 0) {
      p <- p + geom_point(data = singles, position = position_dodge(0.8),
                           shape = 21, size = 4, colour = "black", stroke = 0.8,
                           show.legend = FALSE)
    }
    p
  })

  output$condition_stats <- renderDT({
    stats <- condition_stats_df()
    req(stats)
    datatable(stats, rownames = FALSE,
              options = list(pageLength = 20, dom = "tp")) %>%
      formatStyle("significant", target = "row",
                  backgroundColor = styleEqual(TRUE, "#d4edda"))
  })

  # ── 9. Reproducibility ────────────────────────────────────────────────────

  prov_params <- reactive({
    list(
      k1            = input$k1 %||% 10L,
      k2_min        = input$k2_min %||% 3L,
      k2_max        = input$k2_max %||% 10L,
      loo_mode      = input$loo_mode %||% "count",
      loo_n         = (switch(input$loo_mode %||% "count",
                              count = input$loo_n_count,
                              pct   = input$loo_n_pct,
                              group = input$loo_n_group)) %||% 1,
      agg_fn        = input$agg_fn %||% "median",
      final_k2      = input$final_k2_override %||% rv$optimal_k2 %||% NA,
      optimal_k2    = rv$optimal_k2,
      samples       = rv$sample_names,
      condition_col = rv$condition_col,
      python_backend = PY_AVAILABLE,
      run_date      = as.character(Sys.time()),
      r_version     = paste(R.version$major, R.version$minor, sep = ".")
    )
  })

  output$session_info_out <- renderPrint({
    p <- prov_params()
    cat("── NeighborhoodR Analysis Parameters ──\n\n")
    for (nm in names(p))
      cat(sprintf("  %-20s %s\n", paste0(nm, ":"), paste(p[[nm]], collapse = ", ")))
    cat("\n── Python backend active:", p$python_backend, "──\n")
    if (PY_AVAILABLE) {
      cat("sklearn version:", reticulate::py_eval("__import__('sklearn').__version__"), "\n")
      cat("numpy  version:", reticulate::py_eval("__import__('numpy').__version__"), "\n")
    }
    cat("\n── R Session Info ──\n\n")
    print(sessionInfo())
  })

  replay_script <- reactive({
    p      <- prov_params()
    snames <- p$samples %||% "sample1"
    ct_cols_v <- vapply(snames, function(s) {
      id <- paste0("ct_col_", make.names(s))
      input[[id]] %||% "celltype"
    }, character(1))

    glue::glue(r"(
# NeighborhoodR Replay Script
# Generated: {p$run_date}
# Python backend was active: {p$python_backend}
# ─────────────────────────────────────────────────────────────────────────────

library(SpatialExperiment); library(RANN); library(dplyr)
USE_PYTHON <- {p$python_backend}
if (USE_PYTHON) library(reticulate)

K1       <- {p$k1}
K2_MIN   <- {p$k2_min}
K2_MAX   <- {p$k2_max}
LOO_MODE <- "{p$loo_mode}"
LOO_N    <- {p$loo_n}
AGG_FN   <- "{p$agg_fn}"
FINAL_K2 <- {p$final_k2 %||% "NULL"}
SEED     <- 42L

# ── Load samples — replace paths ──
spe_paths <- c(
  {paste(paste0('"', snames, '" = "/path/to/', snames, '.rds"'), collapse = ",\n  ")}
)
spe_list <- lapply(spe_paths, readRDS)

ct_cols <- c(
  {paste(paste0('"', snames, '" = "', ct_cols_v, '"'), collapse = ",\n  ")}
)

# ── Source Python helpers (copy PY_CODE block from app.R) ──
# if (USE_PYTHON) reticulate::py_run_string(PY_CODE)

# ── Build niche matrix and run sweep ──
# See app.R::build_niche_data() and .r_loo_stability_sweep() for full source.

# ── Final assignment ──
if (USE_PYTHON) {
  km_res      <- reticulate::py$final_kmeans(niche_mat, FINAL_K2)
  assignments <- as.integer(unlist(km_res$labels))
} else {
  set.seed(SEED)
  km          <- kmeans(niche_mat, centers = FINAL_K2, nstart = 25L, iter.max = 300L)
  assignments <- km$cluster
}
ptr <- 1L
for (s in names(spe_list)) {
  n <- ncol(spe_list[[s]])
  colData(spe_list[[s]])$neighborhood <- paste0("N", assignments[ptr:(ptr + n - 1L)])
  colData(spe_list[[s]])$sample        <- rep(s, n)
  # Per-sample cell-type columns aren't guaranteed to share a name across
  # samples (e.g. cluster-label columns embed the per-sample cluster count),
  # so normalise into a common field for cross-sample use.
  colData(spe_list[[s]])$celltype <- as.character(colData(spe_list[[s]])[[ct_cols[[s]]]])
  ptr <- ptr + n
}
# spe_list keeps full marker data — save as-is for downstream analysis.
saveRDS(spe_list, "neighborhoodR_full_list.rds")

# Concatenate into a single SPE for plotting/summary only. cbind() requires
# identical row counts AND identical colData column names across objects.
# Marker panels and original colData schemas can both differ per sample —
# neither of which this app's own logic reads — so keep only the uniform
# fields every sample now has.
spe_light <- lapply(spe_list, function(s) {
  cd <- colData(s)[, intersect(c("sample", "neighborhood", "celltype"), colnames(colData(s))), drop = FALSE]
  SpatialExperiment(
    assays        = list(placeholder = matrix(numeric(0), nrow = 0, ncol = ncol(s))),
    colData       = cd,
    spatialCoords = spatialCoords(s)
  )
})
joint_spe <- do.call(cbind, unname(spe_light))
saveRDS(joint_spe, "neighborhoodR_joint.rds")
message("Done — full list saved to neighborhoodR_full_list.rds, joint SPE saved to neighborhoodR_joint.rds")
)")
  })

  output$prov_preview <- renderPrint({
    cat(jsonlite::toJSON(prov_params(), pretty = TRUE, auto_unbox = TRUE))
  })

  output$download_script <- downloadHandler(
    filename = function() paste0("neighborhoodR_replay_", Sys.Date(), ".R"),
    content  = function(file) writeLines(replay_script(), file)
  )

  output$download_prov_json <- downloadHandler(
    filename = function() paste0("neighborhoodR_provenance_", Sys.Date(), ".json"),
    content  = function(file)
      writeLines(jsonlite::toJSON(prov_params(), pretty = TRUE, auto_unbox = TRUE), file)
  )

  # ── Update condition map without re-running assignment ────────────────────
  observeEvent(input$update_condition, {
    req(rv$joint_spe)
    cond_col <- input$condition_col
    if (!is.null(cond_col) && nchar(cond_col) &&
        cond_col %in% colnames(colData(rv$joint_spe))) {
      rv$condition_col <- cond_col
      rv$condition_map <- NULL
      rv$has_condition <- TRUE
    } else {
      cmap <- setNames(
        vapply(rv$sample_names, function(s) {
          trimws(input[[paste0("cond_map_", make.names(s))]] %||% "")
        }, character(1)),
        rv$sample_names
      )
      if (any(nchar(cmap) > 0)) {
        rv$condition_col <- NULL
        rv$condition_map <- as.list(cmap)
        snames_v  <- as.character(colData(rv$joint_spe)$sample)
        cond_vec  <- rep(NA_character_, ncol(rv$joint_spe))
        for (sname in rv$sample_names) {
          lbl <- cmap[[sname]]
          if (nchar(lbl)) cond_vec[snames_v == sname] <- lbl
        }
        colData(rv$joint_spe)$condition_map <- cond_vec
        # Keep the full-marker export list in sync with the same labels
        if (!is.null(rv$spe_list_full)) {
          for (sname in rv$sample_names) {
            lbl <- cmap[[sname]]
            colData(rv$spe_list_full[[sname]])$condition_map <-
              rep(if (nchar(lbl)) lbl else NA_character_, ncol(rv$spe_list_full[[sname]]))
          }
        }
        rv$has_condition <- TRUE
      } else {
        rv$has_condition <- FALSE
      }
    }
    showNotification("Condition map updated — visualisation plots refreshed.", type = "message")
  })

  # ── Reset / clear memory ───────────────────────────────────────────────────
  observeEvent(input$reset_session, {
    for (nm in c("spe_list", "sample_names", "niche_mat", "joint_spe",
                 "spe_list_full",
                 "cell_types_v", "ct_encoded_v", "sample_labels",
                 "samp_encoded", "all_cts", "sweep_res", "assignments",
                 "condition_map", "mem_sweep")) {
      rv[[nm]] <- NULL
    }
    rv$has_condition <- FALSE
    rv$is_python     <- FALSE
    rv$optimal_k2    <- NULL
    rv$n_samples     <- NULL
    if (PY_AVAILABLE) tryCatch(reticulate::py$deep_clean(), error = function(e) NULL)
    invisible(gc(full = TRUE, reset = TRUE))
    session$reload()
  })

  # Cleanup when browser tab closes
  session$onSessionEnded(function() {
    if (PY_AVAILABLE) tryCatch(reticulate::py$deep_clean(), error = function(e) NULL)
    gc(full = TRUE)
  })
}

shinyApp(ui, server)
