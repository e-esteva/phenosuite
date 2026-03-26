# PhenoSuite

An integrated bioinformatics platform for spatial and single-cell genomics analysis, built on R Shiny. Developed by the Translational Immunology Cennter (TrIC) at NYU, PhenoSuite provides modular tools for multiplex immunofluorescence imaging, cellular phenotyping, spatial transcriptomics, and automated AI-assisted classification.

## Modules

| Module | Description |
|--------|-------------|
| **analysis_engine_multi_config_streamline** | Multi-sample workflow orchestration with flexible configuration |
| **analysis_pre** | Pre-processing tools: segmentation file analysis, spatial gating, and SpatialExperiment analysis |
| **automated_phenotyping** | GPT-powered automated cell phenotyping and classification |
| **circos-artist** | Refine and style circos plot visualizations |
| **circos-builder** | Generate circos plots for cell-cell interactions |
| **interactive_gating** | Draw polygonal gates on 2D scatter plots |
| **masquerade** | Cell-mask overlay builder for multiplex TIFF images |
| **merging_integration** | Merge multiple analysis outputs |
| **modify_clusters** | Rename and reclassify cell clusters |
| **modify_spe** | Modify SpatialExperiment metadata and annotations |
| **pcf-builder** | Proximity Cell Frequency analysis |
| **pcf-v2** | Proximity Cell Frequency analysis (v2) |
| **spatialExploreR** | General-purpose spatial data upload, visualization, and mutation |
| **spatial_interactions** | Model and visualize spatial cell interactions |

## Tech Stack

- **R / Shiny** — primary language and web framework
- **Python** — image processing and ML (via `reticulate`)
- **Key R packages**: SpatialExperiment, SingleCellExperiment, Seurat, ComplexHeatmap, circlize, ggplot2, tidyverse, DT, openai, promises/future
- **Key Python libraries**: tifffile, numpy, scipy, scikit-learn

## Project Structure

```
phenosuite/
├── analysis_engine_multi_config_streamline/
├── analysis_pre/
│   ├── segmentation_file_analysis/
│   ├── spatial_gating/
│   └── spe_analysis/
├── automated_phenotyping/
├── circos-artist/
├── circos-builder/          # dev + production
├── interactive_gating/      # dev + production
├── masquerade/
├── merging_integration/
├── modify_clusters/         # dev + production
├── modify_spe/
├── pcf-builder/
├── pcf-v2/
├── spatialExploreR/
├── spatial_interactions/
└── utils/                   # shared utilities
    ├── RunPhenomenalist-shiny/
    ├── spatial-shiny/
    └── gpts/
```

Modules with both `dev/` and `production/` subdirectories support separate development and deployment workflows. Simpler modules use a single `app.R`; more complex ones split into `server.R` and `ui.R`.

## Shared Utilities

The `utils/` directory contains shared R and Python code sourced by multiple modules:

- **RunPhenomenalist-shiny/** — core Phenomenalist analysis functions and workflow wrappers
- **spatial-shiny/** — circos rendering, spatial interaction modeling, heatmap generation, and Vectra imaging libraries (Python)
- **gpts/** — GPT model configuration for automated phenotyping

## Running Locally in RStudio

1. Clone the repository and open the project in RStudio.
2. Install the required R packages listed in the Tech Stack section above.
3. For modules that use Python (e.g., `masquerade`), ensure Python dependencies are installed and accessible via `reticulate`.
4. To run a module, open its entry file in RStudio:
   - For single-file apps: open `app.R` (e.g., `spatialExploreR/app.R`)
   - For two-file apps: open `server.R` or `ui.R` (e.g., `circos-builder/production/server.R`)
5. Click **Run App** in RStudio, or run from the console:
   ```r
   shiny::runApp("spatialExploreR")
   ```

Shared utilities in `utils/` are sourced at runtime by each module. If running outside the default Shiny Server path, update the `source()` paths in the relevant `server.R`/`app.R` files to point to your local `utils/` directory.

## Data Flow

1. **Input** — CSV, RDS, or SpatialExperiment/Seurat objects from segmentation pipelines
2. **Pre-processing** — Clean and annotate data via `analysis_pre` modules
3. **Analysis** — Clustering, gating, spatial interactions, proximity analysis
4. **AI classification** — Optional GPT-powered automated phenotyping
5. **Visualization** — Scatter plots, heatmaps, circos diagrams
6. **Export** — Download results as ZIP, PNG, or CSV
