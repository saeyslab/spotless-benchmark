# Build MuSiC
FROM csangara/seurat:4.1.0
# If old Seurat version, install remotes and R.utils too
RUN R -e "install.packages(c('BiocManager', 'MCMCpack', 'nnls')); BiocManager::install('Biobase')"
RUN R -e "remotes::install_github(c('renozao/xbioc', 'xuranw/MuSiC'), upgrade='never')"