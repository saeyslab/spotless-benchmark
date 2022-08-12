FROM ubuntu:18.04

RUN apt-get update \
    && apt-get install --no-install-recommends --yes \
        curl zip unzip g++ vim wget ca-certificates  git \
        gcc-8 libssl1.0.0 software-properties-common  \
    && add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get update \
    && apt-get install --yes --only-upgrade libstdc++6 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh
RUN /opt/conda/condabin/conda init bash
COPY extra_code.txt /tmp/
RUN /opt/conda/condabin/conda create -n stride python=3.8
ENV PATH /opt/conda/envs/stride/bin:/opt/conda/bin:$PATH
SHELL ["/bin/bash", "-c"]
RUN source activate stride && \
    git clone https://github.com/csangara/STRIDE.git && cd STRIDE && \
    # Add support for H5AD files
    sed -n 21,29p /tmp/extra_code.txt > tmp && sed -i 83rtmp src/STRIDE/ModelTrain.py && \
    sed -n 9,19p /tmp/extra_code.txt > tmp && sed -i 43rtmp src/STRIDE/ModelTrain.py && \
    sed -i '/import read_10X_h5/ s/$/, read_h5ad/' src/STRIDE/ModelTrain.py && \
    sed -n 1,7p /tmp/extra_code.txt > tmp && sed -i -e 35rtmp src/STRIDE/utility/IO.py && \
    sed -i "10i import scanpy" src/STRIDE/utility/IO.py && \
    # borrowed from https://github.com/DongqingSun96/STRIDE/commit/5b7d845178ddecd6c5671b079fc0e6da40288be4
    sed -i "201i \    celltype_topic_norm_df = celltype_topic_norm_df.fillna(0)" src/STRIDE/Deconvolution.py && \
    sed -i "144i \    celltype_topic_norm_df = celltype_topic_norm_df.fillna(0)" src/STRIDE/Deconvolution.py && \
    pip install -r requirements.txt && python setup.py install

RUN echo "source activate stride" >> ~/.bashrc




