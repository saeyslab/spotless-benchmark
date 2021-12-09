# Build MuSiC
FROM csangara/seurat:3.2.3
RUN R -e "install.packages(c('remotes', 'BiocManager', 'MCMCpack', 'nnls', 'R.utils')); BiocManager::install('Biobase')"
RUN R -e "remotes::install_github(c('renozao/xbioc', 'xuranw/MuSiC'), upgrade='never')"