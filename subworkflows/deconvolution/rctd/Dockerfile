# Build RCTD
FROM csangara/seurat:4.1.0
# RUN R -e "install.packages(c('remotes', 'R.utils'))"
RUN R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-03-17')); remotes::install_github('csangara/RCTD', upgrade='never')"