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
