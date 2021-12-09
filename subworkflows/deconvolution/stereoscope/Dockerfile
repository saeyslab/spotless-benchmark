# Build stereoscope
FROM continuumio/miniconda3:latest
SHELL ["/bin/bash", "-c"]
RUN conda create -n stereoscope python=3.7 numba jupyter scanpy numpy pandas matplotlib scikit-learn umap-learn anndata scipy Pillow libgcc --channel bioconda --channel conda-forge
RUN source activate stereoscope && pip install torch loompy && \
 git clone https://github.com/almaan/stereoscope && cd stereoscope && ./setup.py install
RUN printf "source activate stereoscope\nexport LD_LIBRARY_PATH=/opt/conda/envs/stereoscope/lib" >> ~/.bashrc