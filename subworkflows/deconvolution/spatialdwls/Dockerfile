# Build spatialDWLS
FROM csangara/seurat:4.1.0
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    python3.6 python3-pip python3-setuptools python3-dev \
    && pip3 install pandas python-igraph networkx python-louvain leidenalg scikit-learn
RUN R -e "install.packages('quadprog'); remotes::install_github('RubD/Giotto', upgrade='never')"