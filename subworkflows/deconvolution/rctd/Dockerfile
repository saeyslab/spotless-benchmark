FROM csangara/seurat:5.1.0

# Previous version (will probably not work anymore)
# csangara/RCTD@f61e351c067af7c81ea3c20eceef4e28d46bfcb7

# Normal RCTD
# RUN R -e "options(repos = c(CRAN = 'https://p3m.dev/cran/__linux__/jammy/2024-09-20')); remotes::install_github('csangara/RCTD@0a0861e3d1e16014a20e9b743d0e19d3b42231f3', upgrade='never')"

# Visium HD speed up
RUN R -e "options(repos = c(CRAN = 'https://p3m.dev/cran/__linux__/jammy/2024-09-20')); remotes::install_github('csangara/RCTD', ref = 'HD', upgrade='never')"
