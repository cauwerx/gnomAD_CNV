# Docker image for gnomAD_CNV_GD_dCR_analysis.wdl
# Provides: R (with data.table, GenomicRanges, openxlsx) + htslib (tabix)
#
# Build & push, e.g.:
#   docker build -t us.gcr.io/YOUR_PROJECT/gnomad-gd-dcr:latest .
#   docker push us.gcr.io/YOUR_PROJECT/gnomad-gd-dcr:latest
# Then set that image as the `docker` input when launching the workflow on Terra.

FROM bioconductor/bioconductor_docker:RELEASE_3_18

RUN apt-get update && \
    apt-get install -y --no-install-recommends tabix bcftools && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('data.table', 'openxlsx'), repos='https://cloud.r-project.org')" && \
    R -e "if (!requireNamespace('GenomicRanges', quietly = TRUE)) BiocManager::install('GenomicRanges', update = FALSE, ask = FALSE)"
