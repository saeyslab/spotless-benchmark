FROM nvidia/cuda:10.2-base-ubuntu18.04

# install os packages, also update gcc and libstdc++
RUN apt-get update \
    && apt-get install --no-install-recommends --yes \
        curl \
        unzip \
        g++ \
        vim \
        wget \
        ca-certificates \
        git \
        gcc-8 libssl1.0.0 software-properties-common \
    && add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get update \
    && apt-get install --yes --only-upgrade libstdc++6 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# see http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# install miniconda3 - https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh
RUN /opt/conda/condabin/conda init bash
COPY environment.yml /tmp/
RUN /opt/conda/condabin/conda env create -f /tmp/environment.yml && \
    /opt/conda/condabin/conda clean --all --yes --quiet

ENV PATH /opt/conda/envs/tangram-env/bin:/opt/conda/bin:$PATH
SHELL ["/bin/bash", "-c"]
RUN source activate tangram-env && pip install tangram-sc==1.0.3 squidpy==1.2.2
RUN echo "source activate tangram-sc" >> ~/.bashrc
