FROM --platform=linux/amd64 rocker/shiny:4.4.0

# ── Layer 1: System dependencies ──────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    patch \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    libglpk-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libmagick++-dev \
    cmake \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    libpython3.10 \
    libtirpc-dev \
    && rm -rf /var/lib/apt/lists/*

# ── Layer 2: Python virtual environment ───────────────────
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir \
    "numpy<2" "pandas<2" \
    scipy tifffile imagecodecs scikit-learn matplotlib seaborn \
    "rpy2==3.5.17" \
    anndata

# ── Layer 3: Bioconductor packages ────────────────────────
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
# Install matrixStats from current CRAN (Bioconductor 3.20 needs >= 1.4.1)
RUN R -e "install.packages('matrixStats', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install('SummarizedExperiment', ask=FALSE, update=FALSE, force=TRUE)"
RUN R -e "BiocManager::install('SingleCellExperiment', ask=FALSE, update=FALSE, force=TRUE)"
RUN R -e "BiocManager::install('SpatialExperiment', ask=FALSE, update=FALSE, force=TRUE)"
RUN R -e "BiocManager::install(c('ComplexHeatmap','scater','scran','scuttle','MatrixGenerics'), ask=FALSE, update=FALSE, force=TRUE)"
# Verify critical packages installed
RUN R -e "library(SpatialExperiment); library(SingleCellExperiment); cat('Bioconductor OK\n')"

# ── Layer 4: CRAN packages ───────────────────────────────
RUN install2.r --error --skipinstalled \
    broom \
    circlize \
    colorspace \
    colourpicker \
    cowplot \
    curl \
    data.table \
    digest \
    dplyr \
    DT \
    FNN \
    future \
    ggplot2 \
    ggpubr \
    ggrepel \
    ggrastr \
    ggsci \
    ggsignif \
    glue \
    gridExtra \
    harmony \
    hdf5r \
    htmlwidgets \
    igraph \
    janitor \
    jsonlite \
    ks \
    lubridate \
    mclust \
    openai \
    patchwork \
    pheatmap \
    plotly \
    promises \
    RColorBrewer \
    readr \
    remotes \
    reticulate \
    RSpectra \
    scales \
    scattermore \
    Seurat \
    shinyFiles \
    shinyjs \
    shinyWidgets \
    shinycssloaders \
    shinydashboard \
    sp \
    stringi \
    stringr \
    tidyr \
    tidyverse \
    uwot \
    viridis \
    zip

# ── Layer 5b (intentionally omitted): HF/torch backend ───
# The gemma_phenotyper GUI is Ollama-only: it calls a running Ollama server
# over HTTP from R and needs no Python ML stack at all. torch + transformers
# would add ~2-3 GB to this image, and CPU-only inference on a 12B checkpoint
# is impractical anyway (hours per sample).
#
# The HuggingFace backend lives in phenosuite-CLI/gemma-phenotyper, which runs
# on GPU nodes where it belongs. If you ever want it in this container, install:
#   torch (--index-url https://download.pytorch.org/whl/cpu for CPU-only)
#   transformers>=5.10.0   # Gemma 4 floor; >=4.50.0 suffices for Gemma 3 only
#   peft accelerate bitsandbytes

# ── Layer 5: GitHub packages ─────────────────────────────
RUN R -e "remotes::install_github('igordot/phenomenalist', dependencies=FALSE)"

# ── Layer 6: Configuration ───────────────────────────────
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY docker-entrypoint.sh /docker-entrypoint.sh
RUN chmod +x /docker-entrypoint.sh

# ── Layer 7: Provenance metadata ─────────────────────────
ARG GIT_SHA=unknown
ARG IMAGE_DIGEST=unknown
ENV PHENOSUITE_GIT_SHA=${GIT_SHA}
ENV PHENOSUITE_IMAGE_DIGEST=${IMAGE_DIGEST}

EXPOSE 3838

ENTRYPOINT ["/docker-entrypoint.sh"]
