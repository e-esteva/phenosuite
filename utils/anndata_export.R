# AnnData (.h5ad) writer for SpatialExperiment / SingleCellExperiment objects.
#
# Replaces the earlier loom writer. Loom is a frozen format that scanpy/squidpy
# read only via a compatibility shim, and it has no first-class home for spatial
# coordinates — they had to ride along as ordinary column attributes. AnnData is
# what the Python spatial-omics stack (scanpy, squidpy) actually operates on,
# and it stores coordinates properly in obsm["spatial"], which is exactly where
# squidpy looks for them.
#
# Written through the reference Python implementation via reticulate rather than
# by hand: .h5ad is HDF5 with per-element "encoding-type"/"encoding-version"
# attributes that readers validate, and hand-rolling those is a good way to
# produce a file that opens today and breaks on the next anndata release.
#
# Requires: reticulate, and the `anndata` Python package in the interpreter
# reticulate binds to (RETICULATE_PYTHON, set to /opt/venv/bin/python in
# docker-compose.yml).

# spe        : SpatialExperiment or SingleCellExperiment
# file_path  : output .h5ad path
# assay      : assay to store as adata.X. Default picks exprs > logcounts >
#              counts > first available, matching how the apps read intensities.
# obs_extra  : optional named list of per-cell vectors to add to adata.obs on
#              top of colData (e.g. annotations computed after the object was
#              built).
# layers     : if TRUE, every other assay is written to adata.layers so nothing
#              is silently dropped.
write_h5ad <- function(spe,
                       file_path,
                       assay      = NULL,
                       obs_extra  = NULL,
                       layers     = TRUE) {
  require(reticulate)
  require(SummarizedExperiment)

  ad <- tryCatch(reticulate::import("anndata"), error = function(e) NULL)
  if (is.null(ad))
    stop("write_h5ad: the Python 'anndata' package is not importable from ",
         "reticulate's interpreter (", reticulate::py_config()$python, "). ",
         "Install it with: pip install anndata")

  avail <- assayNames(spe)
  if (length(avail) == 0) stop("write_h5ad: object has no assays")
  if (is.null(assay)) {
    assay <- c("exprs", "logcounts", "counts")[
      c("exprs", "logcounts", "counts") %in% avail][1]
    if (is.na(assay)) assay <- avail[1]
  }
  if (!assay %in% avail)
    stop("write_h5ad: assay '", assay, "' not found. Available: ",
         paste(avail, collapse = ", "))

  # AnnData is observations x variables — cells as ROWS. SPE assays are the
  # transpose (features x cells), so every matrix has to be flipped on the way
  # out. Getting this backwards produces a file that loads but is nonsense.
  to_dense_t <- function(m) t(as.matrix(m))

  X <- to_dense_t(assay(spe, assay))

  obs <- as.data.frame(colData(spe))
  if (!is.null(obs_extra)) {
    for (nm in names(obs_extra)) {
      v <- obs_extra[[nm]]
      if (length(v) == ncol(spe)) obs[[nm]] <- v
      else message("write_h5ad: dropping obs_extra['", nm, "'] (length mismatch)")
    }
  }
  # h5ad cannot store list/S4 columns; coerce anything exotic to character so a
  # single awkward column doesn't fail the whole export.
  for (nm in names(obs)) {
    if (!is.atomic(obs[[nm]])) obs[[nm]] <- as.character(obs[[nm]])
    if (is.factor(obs[[nm]]))  obs[[nm]] <- as.character(obs[[nm]])
  }

  var <- as.data.frame(rowData(spe))
  if (ncol(var) == 0) var <- data.frame(row.names = rownames(spe))
  var$marker <- rownames(spe)

  cell_ids <- colnames(spe)
  if (is.null(cell_ids)) cell_ids <- paste0("cell_", seq_len(ncol(spe)))
  rownames(obs) <- cell_ids
  rownames(X)   <- cell_ids
  colnames(X)   <- rownames(spe)
  rownames(var) <- rownames(spe)

  # obsm / layers are built up FIRST and handed to the constructor. Assigning
  # into adata$obsm[[...]] after construction looks like it works from R but
  # does not mutate the underlying Python dict, so the keys silently vanish.
  obsm_list <- list()

  # Spatial coordinates -> obsm["spatial"], the key squidpy expects.
  if (inherits(spe, "SpatialExperiment")) {
    xy <- tryCatch(as.matrix(spatialCoords(spe)), error = function(e) NULL)
    if (!is.null(xy) && nrow(xy) == ncol(spe)) {
      storage.mode(xy) <- "double"
      dimnames(xy) <- NULL
      obsm_list[["spatial"]] <- xy
    }
  }

  # Reduced dimensions (UMAP etc.) -> obsm, using scanpy's X_ prefix convention.
  if (methods::is(spe, "SingleCellExperiment")) {
    for (d in SingleCellExperiment::reducedDimNames(spe)) {
      m <- tryCatch(as.matrix(SingleCellExperiment::reducedDim(spe, d)),
                    error = function(e) NULL)
      if (!is.null(m) && nrow(m) == ncol(spe)) {
        storage.mode(m) <- "double"
        dimnames(m) <- NULL
        obsm_list[[paste0("X_", tolower(d))]] <- m
      }
    }
  }

  layers_list <- list()
  if (isTRUE(layers)) {
    for (a in setdiff(avail, assay)) {
      m <- tryCatch(to_dense_t(assay(spe, a)), error = function(e) NULL)
      if (!is.null(m)) {
        dimnames(m) <- NULL
        layers_list[[a]] <- m
      }
    }
  }

  args <- list(X = X, obs = obs, var = var)
  if (length(obsm_list)   > 0) args$obsm   <- obsm_list
  if (length(layers_list) > 0) args$layers <- layers_list
  adata <- do.call(ad$AnnData, args)

  adata$write_h5ad(file_path)
  message("write_h5ad: wrote ", file_path, "  (", ncol(spe), " cells x ",
          nrow(spe), " markers, X = '", assay, "')")
  invisible(file_path)
}
