# Build Seurat-Disk with Seurat 4.0.5
FROM csangara/seurat:4.0.5
RUN apt-get update -y && apt-get install -y libglpk-dev
RUN R -e "install.packages(c('remotes', 'R.utils', 'hdf5r'))"
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"