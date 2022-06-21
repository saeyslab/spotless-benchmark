# Build MuSiC
FROM csangara/seurat:4.1.0
# If old Seurat version, install remotes and R.utils too
RUN R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-03-17')); install.packages(c('BiocManager', 'MCMCpack', 'nnls')); BiocManager::install('Biobase')"
RUN R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-03-17')); remotes::install_github(c('renozao/xbioc@0.1.19', 'xuranw/MuSiC@7c5834830223957f5d8134c86d6acf653bfff4e7'), upgrade='never')"