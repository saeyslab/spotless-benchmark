# base image maintained by the NVIDIA CUDA Installer Team - https://hub.docker.com/r/nvidia/cuda/
FROM nvidia/cuda:10.2-cudnn7-devel-ubuntu18.04

# NVIDIA changed their GPG key in April 2022: https://developer.nvidia.com/blog/updating-the-cuda-linux-gpg-repository-key/
# Fix by aalphabio: https://gitlab.com/nvidia/container-images/cuda/-/issues/158
RUN apt-key del "7fa2af80" \
&& export this_distro="$(cat /etc/os-release | grep '^ID=' | awk -F'=' '{print $2}')" \
&& export this_version="$(cat /etc/os-release | grep '^VERSION_ID=' | awk -F'=' '{print $2}' | sed 's/[^0-9]*//g')" \
&& apt-key adv --fetch-keys "https://developer.download.nvidia.com/compute/cuda/repos/${this_distro}${this_version}/x86_64/3bf863cc.pub" \
&& apt-key adv --fetch-keys "https://developer.download.nvidia.com/compute/machine-learning/repos/${this_distro}${this_version}/x86_64/7fa2af80.pub"

# install os packages
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

# create conda environment yaml file
COPY environment.yml /tmp/
RUN /opt/conda/condabin/conda config --set channel_priority strict \
    && /opt/conda/condabin/conda env create -f /tmp/environment.yml \
    && echo "source activate cell2loc_env" >> ~/.bashrc \
    && /opt/conda/condabin/conda clean --all --yes --quiet
ENV PATH /opt/conda/envs/cell2loc_env/bin:/opt/conda/bin:$PATH

# install cell2location 
RUN /bin/bash -c "pip install git+https://github.com/BayraktarLab/cell2location.git@42fe3ccda054b09cf737ce038eb00b2cec9ea10a#egg=cell2location[tutorials]"