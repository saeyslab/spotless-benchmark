# Build RCTD
FROM csangara/seurat:3.2.3
RUN R -e "install.packages(c('remotes', 'R.utils'))"
RUN R -e "remotes::install_github('dmcable/RCTD', upgrade='never')"