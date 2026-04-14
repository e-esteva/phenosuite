# PhenoSuite

An integrated bioinformatics platform for **single-cell spatial omics** analysis, built on R Shiny. Developed by the Translational Immunology Center (TrIC) at NYU, PhenoSuite provides modular tools for multiplexed immunofluorescence imaging (CODEX/PhenoCycler), spatial transcriptomics (MERFISH), cellular phenotyping, spatial interaction analysis, and multi-modal data integration.

The entire platform ships as a **single Docker image** — clone the repo, run one command, and browse to `localhost:3838`. No institutional server access needed.

---

## Quick Start (Docker)

### Prerequisites

1. **Docker Desktop** — download and install from [docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop/)
2. **Git** — to clone the repository
3. **~10 GB free RAM** allocated to Docker (Settings → Resources → Memory)
4. **~8 GB free disk space** for the Docker image

### Setup

```bash
# 1. Clone the repository
git clone https://github.com/e-esteva/phenosuite.git
cd phenosuite

# 2. Build and start the container (first build: 30–60 min)
docker compose up --build

# 3. Open your browser to:
#    http://localhost:3838/site/
```

That's it. The landing page will load with all apps accessible.

> **First build is slow.** Compiling Bioconductor packages (SpatialExperiment, ComplexHeatmap, etc.) and Seurat takes 30–60 minutes. Subsequent builds use Docker's layer cache and are nearly instant.

> **Using `docker-compose` instead of `docker compose`?** If you have the legacy `docker-compose` binary (hyphenated), it works the same way: `docker-compose up --build`.


### Environment variables

Optional — copy `.env.example` to `.env` to set:

```bash
OPENAI_API_KEY=sk-...   # required only for the Automated Phenotyping (GPT) app
```

---

## Platform Overview

PhenoSuite is organized into **four analysis modules** exposed via a tabbed landing page:

### Analysis Engine
End-to-end workflow for loading and pre-processing spatial data.

| App | Description |
|---|---|
| **Launch Analysis Engine** | Multi-sample workflow orchestration |
| **Spatial Gating of Segmentation File** | Draw spatial gates on raw segmentation CSVs |
| **Compare Segmentation Files** | Side-by-side QC of multiple segmentation outputs |
| **MERFISH Spatial Transcriptomics** | Full MERFISH pipeline: QC → normalize → cluster → DE → SVG → export SPE |
| **Compare SPE Objects** | Harmonized comparison across SpatialExperiment files |
| **Merge / Integrate SPE Objects** | Concatenate or Harmony-integrate multiple SPEs |
| **Explore SPE Objects** | Interactive browser for any SpatialExperiment |
| **CODEX + MERFISH Integration** | Multi-modal spatial integration (see below) |

### Phenotyping Module
Cluster annotation and cell-type refinement.

| App | Description |
|---|---|
| **Sub-cluster / Re-cluster** | Re-run clustering on subsets of an existing SPE |
| **Interactive Gating** | Polygonal 2D gating on scatter plots |
| **Automated Phenotyping** | GPT-powered LLM-driven cluster annotation |
| **Manually Annotate Clusters** | Live heatmap-based manual cluster labeling |
| **Masquerade** | Generate cluster-overlay masks on multiplexed TIFFs |

### Spatial Module
Cellular neighborhood and interaction analysis.

| App | Description |
|---|---|
| **Pair Correlation Function (PCF)** | PCF analysis across cell-type pairs |
| **Pairwise Log-Odds Interactions** | Enrichment/depletion of cell-type co-occurrence |

### Graphics Module
Publication-ready visualization.

| App | Description |
|---|---|
| **PCF Builder** | Multi-sample PCF comparison with shared reference |
| **Circos Artist** | Single-sample circos plot styling |
| **Circos Builder** | Multi-sample circos harmonization |

---

## Automated Phenotyping — Prompt Algorithms

The **Automated Phenotyping** app uses a large language model (GPT) to assign cell-type labels to each cluster in a `SpatialExperiment`. For every cluster, the app computes per-marker mean expression across its member cells, picks the most informative markers, and asks the LLM to name the cell type. The exact way markers are selected and the prompt is phrased is controlled by the **Prompt algorithm** dropdown. Two algorithms are currently available, and both can be optionally grounded in a tissue context via the **Tissue type** text box (e.g. `spleen`, `tonsil`, `lung`). The text box is always shown — if left blank, no tissue clause is added to the prompt.

### v1 — Symmetric single choice

**Marker selection.** For each cluster, the app takes markers with mean expression above the 95th percentile *and* below the 5th percentile of that cluster's marker distribution. The thresholds are symmetric around the median (top 5% and bottom 5%), hence *symmetric*.

**Encoding.** Markers are passed to the model as raw `marker:value` pairs — the underlying numeric expression is preserved, so the model sees quantitative evidence rather than a discretised `+/-` call.

**Prompt shape.** The model is asked a single question and instructed to return a three-word label. No alternatives, no explanation:

> *"what celltype is described by {markers}? [This is a {tissue}.] Give me a 3 word response"*

**Output.** Exactly one cell-type label per cluster. Fast, low-token, easy to parse. Good when markers are unambiguous and you trust the model's first call.

### v2 — Asymmetric 2-choice

**Marker selection.** Markers above the 90th percentile (top 10%, a wider "high" set) and below the 5th percentile (bottom 5%) — asymmetric thresholds that surface more positive evidence than negative.

**Encoding.** Markers are passed as `marker+` / `marker-` tokens (thresholded at zero), mimicking the `CD4+CD8-` shorthand immunologists use. The model sees a discretised signature rather than raw values.

**Prompt shape.** The model is asked to propose **two candidate cell types** and justify each one, forced into a structured format so the app can parse both choices:

> *"what celltype is described by {markers}? [This is a {tissue}.] Give me a 3 word response. Give two choices. Explain. Format response like: Choice X: choice; Explanation X: explanation."*

**Output.** Two candidate labels plus a free-text explanation for each, stored alongside the cluster. Both labels are retained in the annotation (joined with a comma), and the explanations are preserved in the downloaded report for auditing. Use this when markers are ambiguous, when you want a second opinion per cluster, or when you plan to manually pick between the two candidates afterwards.

### Post-processing (both algorithms)

After per-cluster calls, the app runs two additional LLM passes that are independent of the chosen prompt algorithm:

1. **Harmonisation.** Near-duplicate labels (`CD4 T cell` vs `CD4T cells`, `NK cell` vs `Natural Killer cell`) are collapsed to a canonical form via a single JSON remapping call. Title Case is enforced; biologically distinct populations are kept separate.
2. **Broad lineage.** Harmonised labels are mapped to coarse immune lineage buckets (`T Cell`, `B Cell`, `NK Cell`, `Macrophage`, `Dendritic Cell`, `Neutrophil`, `Mast Cell`, `Monocyte`) for lineage-level plots. Non-immune labels (stromal, endothelial, epithelial) are preserved as-is.

Both steps are automatic and run after the per-cluster annotation finishes — the choice of `v1` vs `v2` only affects the per-cluster labeling call.

### Choosing between them

| | v1 (symmetric) | v2 (asymmetric 2-choice) |
|---|---|---|
| Markers used | top 5% / bottom 5% | top 10% / bottom 5% |
| Marker encoding | raw `marker:value` | thresholded `marker+` / `marker-` |
| Answers per cluster | 1 | 2 + explanations |
| Token cost per call | lower | ~2–3× higher |
| Best when | panels are clean, labels unambiguous | panels are large, markers overlap, or you want a second opinion |

---

## Provenance Sidecar (Reproducibility)

Every PhenoSuite app bundles a **provenance sidecar** into its output ZIP. The sidecar is a pair of files — `provenance.json` and `replay.R` — written to the session temp dir alongside the normal analysis outputs, so they travel with the download automatically. Together they capture enough information to audit, re-run, or cite an analysis after the fact without needing access to the original Shiny session.

### What's captured

`provenance.json` records:

| Field | Contents |
|---|---|
| `app_name` | Which PhenoSuite app ran (`masquerade`, `merfish`, `automated_phenotyping`, …) |
| `session_id` | Shiny session token — unique per browser tab |
| `timestamps` | `session_start`, `analysis_start`, `analysis_end` (ISO-8601 with timezone) |
| `environment.r_version` | Full R version string from `sessionInfo()` |
| `environment.platform` | OS / architecture |
| `environment.docker_image_digest` | Docker image digest (from `PHENOSUITE_IMAGE_DIGEST` build arg) |
| `environment.git_sha` | Git commit SHA of the repo at build time (from `PHENOSUITE_GIT_SHA`) |
| `environment.packages` | Every attached and loaded R package, with version |
| `inputs[]` | For each uploaded file: original name, size in bytes, **SHA-256** hash |
| `parameters` | Every Shiny input control (sliders, dropdowns, text boxes, checkboxes) captured at analysis-start time, with `fileInput` entries stripped |
| `seeds` | Named random seeds used during the run (e.g. `clustering_leiden`, `downsampling`) |
| `outputs[]` | Every file written to the temp dir, with relative path, size, and SHA-256 |
| `custom_metadata` | App-specific extras — e.g. `interactive_gating` stores the full gate tree, `analysis_engine_multi_config` stores the per-sample config table, `spatial_gating` stores the drawn polygon geometry |

`replay.R` is an auto-generated, self-contained R script that:

1. `library()`-loads every package that was attached during the original run
2. Re-seeds any random generators via the `seeds` table
3. Declares a `params <- list(...)` block pre-populated with the captured UI values
4. Declares an `input_files <- c(...)` block with the original file names
5. **Verifies every input file's SHA-256** against the hash in the sidecar before proceeding — mismatched inputs produce a warning with the expected/actual digests so you know immediately if a file has been modified or corrupted
6. Creates a `replay_output/` directory and, for apps with a registered replay template, calls the relevant `phenomenalist` / `RunPhenomenalist-shiny` function with the captured parameters. Apps without a template get a commented block listing the params and input files so a human can finish the script

The script is written to be runnable via `Rscript replay.R` inside any PhenoSuite container — the Docker image digest in the sidecar tells you exactly which image version to pull.

### Which apps produce it

The sidecar is enforced across **every app served from the landing page**. The table below lists each app, which touch points are wired, and any app-specific notes:

| Module | App | Inputs hashed | Params captured | Outputs hashed | Notes |
|---|---|---|---|---|---|
| **Analysis Engine** | Launch Analysis Engine (`analysis_engine_multi_config_streamline`) | ✓ (all uploaded segmentation / config files) | ✓ | ✓ | Per-sample config table stored in `custom_metadata.sample_configs`; known seeds pre-recorded |
| | Spatial Gating of Segmentation File (`analysis_pre_processing/spatial_gating`) | ✓ | ✓ | ✓ | Drawn polygon vertices stored in `custom_metadata.gates`; written on **Download Gated Cells** |
| | Compare Segmentation Files (`analysis_pre_processing/segmentation_file_analysis`) | ✓ (multi-file) | ✓ | ✓ | Plots are written directly into the provenance dir; download handler zips the full dir |
| | Compare SPE Objects (`analysis_pre_processing/spe_analysis`) | ✓ (multi-RDS) | ✓ | ✓ | Same pattern as segmentation_file_analysis |
| | MERFISH Spatial Transcriptomics (`merfish`) | ✓ (expression + metadata) | ✓ | ✓ | Dedicated `dl_provenance` download button because outputs are many individual files |
| | Merge / Integrate SPE Objects (`merging_integration`) | ✓ (all SPEs) | ✓ | ✓ | Sidecar copied into the merged RDS bundle |
| | Explore SPE Objects (`spatialExploreR`) | ✓ | ✓ | ✓ | Dedicated `dl_provenance` download button |
| | CODEX + MERFISH Integration (`multimodal_integration`) | ✓ (CODEX + MERFISH SPEs) | ✓ | ✓ | Dedicated `dl_provenance` download button |
| **Phenotyping** | Sub-cluster / Re-cluster (`modify_spe`) | ✓ | ✓ | ✓ | Two independent trackers — one per sub-clustering / re-clustering panel; known seeds pre-recorded |
| | Sub-clustering (legacy `sub_clustering`) | ✓ | ✓ | ✓ | Still enforced even though not linked from the homepage |
| | Interactive Gating (`interactive_gating/production`) | ✓ | ✓ (captured at gate-apply time) | ✓ | Full gate tree (hierarchy + thresholds) stored in `custom_metadata.gate_tree`; sidecar finalised inside the **Download All** handler |
| | Automated Phenotyping (`automated_phenotyping`) | ✓ | ✓ | ✓ | Captures GPT prompt algorithm, model selection, tissue context, and per-cluster usage reports |
| | Manually Annotate Clusters (`modify_clusters/dev`) | ✓ | ✓ | ✓ | Dev variant is the one linked from the homepage; production variant is also wired |
| | Masquerade (`masquerade`) | ✓ (image + spatial metadata + optional marker whitelist) | ✓ | ✓ | Sidecar bundled into the download ZIP alongside the TIFF |
| **Spatial** | Pair Correlation Function (`pcf-v2`) | ✓ (multi-CSV) | ✓ | ✓ | |
| | Pairwise Log-Odds Interactions (`spatial_interactions`) | ✓ | ✓ | ✓ | |
| **Graphics** | PCF Builder (`pcf-builder`) | ✓ | ✓ | ✓ | |
| | Circos Artist (`circos-artist`) | ✓ | ✓ | ✓ | |
| | Circos Builder (`circos-builder/dev`) | ✓ (multi log-odds CSV) | ✓ | ✓ | Dev variant is linked from the homepage; production variant is also wired |

Every app uses the same `ProvenanceTracker` R5 class defined in `utils/provenance.R`. Touch points are deliberately small (source + init + `register_input` + `capture_parameters` + `analysis_started` + `analysis_completed`), so adding the sidecar to a new app is ~10 lines of boilerplate.

### Using the sidecar

```r
# Inside any replay environment:
library(jsonlite)
prov <- fromJSON("provenance.json", simplifyVector = FALSE)

# Check which image produced this output
prov$environment$docker_image_digest
prov$environment$git_sha

# Verify outputs haven't been tampered with
for (out in prov$outputs) {
  actual <- digest::digest(file = out$relative_path, algo = "sha256")
  stopifnot(actual == out$sha256)
}

# Re-run the analysis headlessly
Rscript replay.R
```

### Configuring the image digest / git SHA

To make replay bulletproof, pass the image digest and git SHA at build time so they end up baked into every sidecar:

```bash
docker build \
  --build-arg GIT_SHA=$(git rev-parse HEAD) \
  --build-arg IMAGE_DIGEST=$(docker image inspect phenosuite:latest --format '{{index .RepoDigests 0}}') \
  -t phenosuite:latest .
```

If the build args are unset the sidecar still works — the corresponding fields just read `"unknown"`.

---

## Multi-Modal Integration (CODEX + MERFISH)

A dedicated Shiny app integrates CODEX protein and MERFISH transcript data from the same tissue:

1. **Import & Align** — Upload both SPE objects; apply translation, scaling, and Y-flip to register MERFISH coordinates to CODEX space.
2. **Cell Matching** — Nearest-neighbor spatial matching with configurable distance threshold. Supports CODEX → MERFISH, MERFISH → CODEX, or mutual nearest neighbors.
3. **Integration** — Weighted feature concatenation → joint UMAP → Leiden/Louvain clustering on the shared embedding.
4. **Export** — Download an integrated `SpatialExperiment` object (CODEX in main assay, MERFISH in `altExp("merfish")`) compatible with all downstream Phenomenalist modules.

---

## Tech Stack

- **R 4.4** + **Shiny Server** (via `rocker/shiny` base image)
- **Python 3** virtualenv (for `reticulate`-based image processing)
- **Bioconductor 3.20** — SpatialExperiment, SingleCellExperiment, SummarizedExperiment, ComplexHeatmap, scater, scran, scuttle
- **CRAN** — Seurat, ggplot2, plotly, DT, tidyverse, uwot, igraph, RANN, future, promises, openai, reticulate, and ~50 others
- **GitHub** — [`igordot/phenomenalist`](https://github.com/igordot/phenomenalist) R package
- **Container runtime** — Docker + Docker Compose
- **Deployment** — Shiny Server Open Source v1.5

Full package list: see `Dockerfile`.

---

## Project Structure

```
phenosuite/
├── Dockerfile                       # R + Bioconductor + Python build
├── docker-compose.yml               # Single-service orchestration
├── docker-entrypoint.sh             # Creates symlinks for hardcoded paths
├── shiny-server.conf                # Shiny Server routing
├── .env.example                     # Environment variable template
│
├── site/                            # Static landing page
│   ├── index.html
│   └── assets/
│
├── analysis_engine_multi_config_streamline/
├── analysis_pre_processing/
│   ├── segmentation_file_analysis/
│   ├── spatial_gating/
│   └── spe_analysis/
├── automated_phenotyping/
├── circos-artist/
├── circos-builder/                  # dev + production
├── interactive_gating/              # dev + production
├── masquerade/
├── merfish/                         # MERFISH spatial transcriptomics
├── merging_integration/
├── modify_clusters/                 # dev + production
├── modify_spe/
├── multimodal_integration/          # CODEX + MERFISH integration
├── pcf-builder/
├── pcf-v2/
├── spatialExploreR/
├── spatial_interactions/
├── sub_clustering/
│
└── utils/                           # shared utilities
    ├── RunPhenomenalist-shiny/
    ├── spatial-shiny/
    └── gpts/
```

Apps with both `dev/` and `production/` subdirectories have separate development and release versions. Simpler apps use a single `app.R`; more complex ones split into `server.R` and `ui.R`.

---

## How Docker Fits the Repo

The container uses a single-image approach:

- The repo is **bind-mounted** at `/srv/shiny-server/phenosuite` — any code edit on the host shows up immediately after a `docker compose restart` (no rebuild).
- The entrypoint creates **symlinks** so hardcoded `source('/srv/shiny-server/phenomenalist/...')` and `source('/srv/shiny-server/Phenoptics-Menu/...')` calls in legacy R code resolve correctly without modification.
- **Shiny Server** auto-discovers every subdirectory containing an `app.R` or `server.R`/`ui.R` and serves it as a Shiny app.
- The static landing page is served from `site/` — navigate to `http://localhost:3838/site/` to start.

---

## Data Flow

1. **Input** — Segmentation CSVs, RDS files with SpatialExperiment/Seurat objects, or raw MERFISH expression matrices
2. **Pre-processing** — QC, normalization, spatial gating
3. **Clustering & annotation** — Manual gating, automated (GPT), or clustering-based
4. **Spatial analysis** — PCF, log-odds interactions, circos visualization
5. **Multi-modal integration** — Optional CODEX + MERFISH joint analysis
6. **Export** — RDS (SpatialExperiment), CSV metadata, PDF figures

---
### Common commands

```bash
# Start in background (detached mode)
docker compose up -d

# View logs
docker compose logs -f phenosuite

# Stop the container
docker compose down

# Restart (picks up code changes — repo is bind-mounted)
docker compose restart

# Shell into the running container
docker compose exec phenosuite bash

# View Shiny app-level logs (useful for debugging a crashed app)
docker compose exec phenosuite bash -c "cat /var/log/shiny-server/*.log"
```
## Troubleshooting

**Build fails with "unknown command: docker compose"**
Docker Desktop isn't running, or the Compose plugin isn't installed. Open Docker Desktop from Applications and wait for the whale icon in the menu bar to stop animating, then retry.

**App loads the landing page but individual apps error out**
Check the Shiny Server app-level logs:
```bash
docker compose exec phenosuite bash -c "cat /var/log/shiny-server/*.log"
```
This shows R package errors, missing source files, and runtime exceptions.

**"Terminated" messages in the logs / apps freezing on large datasets**
Shiny worker processes are being OOM-killed. Raise Docker Desktop's memory allocation (Settings → Resources → Memory) to at least 10 GB. The `docker-compose.yml` reserves up to 8 GB for the container.

**Apple Silicon (M1/M2/M3)**
Most packages compile natively on ARM. If a specific package fails, add `platform: linux/amd64` under the `phenosuite:` service in `docker-compose.yml` to fall back to x86 emulation (slower but universally compatible).

**Changed code doesn't show up**
Code changes are picked up automatically because the repo is bind-mounted, but Shiny needs to restart its worker. Either reload the app page in your browser (new workers spawn per session) or run `docker compose restart`.

---

## Docker Usage (Advanced)

### Routing session temp storage to a fast, high-capacity disk

Several PhenoSuite apps write large intermediate files during a session — multi-gigabyte TIFF masks from **Masquerade**, cached spatial experiment objects, PCF/log-odds outputs, etc. By default these live on Docker's root filesystem, which may be small on laptops or shared servers. If you have a dedicated SSD (or any host path with plenty of free space), you can bind-mount it into the container and point every app's temp directory at it.

Edit `docker-compose.yml`:

```yaml
services:
  phenosuite:
    build:
      context: .
      dockerfile: Dockerfile
    ports:
      - "3838:3838"
    volumes:
      - .:/srv/shiny-server/phenosuite
      # Replace the left-hand path with any directory on a fast, high-capacity disk.
      # On Linux, a dedicated SSD mount works well (e.g. /mnt/fast-ssd/phenosuite-sessions).
      # On macOS/Windows, any host path accessible to Docker Desktop is fine.
      - /mnt/fast-ssd/phenosuite-sessions:/apps/home/rtmp
    environment:
      - OPENAI_API_KEY=${OPENAI_API_KEY:-}
      - RETICULATE_PYTHON=/opt/venv/bin/python
      - R_MAX_VSIZE=8Gb
      # PhenoSuite-aware apps (masquerade, modify_spe, automated_phenotyping, …)
      # read this variable to place per-session temp dirs on the mounted SSD.
      - PHENOSUITE_TMPDIR=/apps/home/rtmp
      # TMPDIR redirects R's built-in tempdir() for every other app that uses
      # the default (tempfile(), write.csv to tempdir(), etc.).
      - TMPDIR=/apps/home/rtmp
    shm_size: "4g"
    mem_limit: "8g"
    memswap_limit: "12g"
```

**How it works:**

- `PHENOSUITE_TMPDIR` is read by apps that have been updated to support it. They fall back to R's `tempdir()` when the variable is unset, so the configuration also works out of the box without any bind-mount.
- `TMPDIR` is the standard POSIX/R temp directory variable — setting it reroutes every other app (including anything spawned by `reticulate`) to the same disk.
- Both variables point at `/apps/home/rtmp` inside the container, which is the mount point for your host SSD.

**Host path requirements:**

- The host directory must exist before `docker compose up` — Docker will bind-mount it as-is.
- The `shiny` user inside the container runs as **UID 999**. Make the host directory writable by that UID:
  ```bash
  sudo mkdir -p /mnt/fast-ssd/phenosuite-sessions
  sudo chown 999:999 /mnt/fast-ssd/phenosuite-sessions
  ```
  On macOS/Windows with Docker Desktop, the file-sharing layer handles permissions automatically — just make sure the path is listed under Settings → Resources → File Sharing.

**When this matters most:**

- **Masquerade** — generates full-resolution cluster masks on multiplexed TIFFs; a single session can produce 2–10 GB of output.
- **MERFISH / Spatial Exploration / Multi-modal integration** — large `SpatialExperiment` RDS files are written and re-read during a run.
- **Long-running analyses** — intermediate state accumulates across steps; running out of temp space mid-run aborts the job with a cryptic error.

If you stick with the default configuration (no bind-mount), a Docker-managed named volume (`session_tmp`) is used instead. That volume lives on Docker's own storage and is fine for small-to-medium datasets.

### Giving Shiny more memory

The `mem_limit` / `memswap_limit` / `shm_size` keys in `docker-compose.yml` cap how much RAM the container can use. If Shiny workers keep getting OOM-killed on large datasets (look for `Terminated` in the logs), raise these limits and make sure Docker Desktop itself has been granted at least that much memory under Settings → Resources.

---

## Running Locally Without Docker (Advanced)

If you prefer running individual apps in RStudio:

1. Install R 4.4 and the packages listed in the `Dockerfile`
2. Install Bioconductor 3.20 packages: `SpatialExperiment`, `SingleCellExperiment`, `ComplexHeatmap`, `scater`, `scran`, `scuttle`
3. Install the phenomenalist package: `remotes::install_github("igordot/phenomenalist")`
4. Create symlinks (or update `source()` calls) so `/srv/shiny-server/phenomenalist/` resolves to your repo root
5. Open any `app.R` or `server.R` in RStudio and click **Run App**

The Docker setup is strongly recommended — it handles all of the above automatically.

---

## License & Citation

See the Phenomenalist package and individual module headers for license information. Please cite the Translational Immunology Center (TrIC) at NYU when using this platform in publications.
