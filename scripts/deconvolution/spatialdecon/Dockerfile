# Build RCTD
FROM csangara/seurat:r4.1
RUN R -e "install.packages(c('remotes', 'BiocManager'))"
RUN R -e "remotes::install_github(c('Nanostring-Biostats/GeomxTools', 'Nanostring-Biostats/SpatialDecon'), \
                build_vignettes = FALSE, upgrade = 'never')"