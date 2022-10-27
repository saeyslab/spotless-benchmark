FROM csangara/seurat:4.1.0
#ENV PATH="/root/miniconda3/bin:${PATH}"
#ARG PATH="/root/miniconda3/bin:${PATH}"
RUN R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-06-17')); install.packages(c('philentropy', 'BiocManager')); BiocManager::install('scuttle')"
RUN apt-get update && apt-get install -y wget git && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

RUN /opt/conda/condabin/conda init bash
RUN /opt/conda/condabin/conda create -n dstg python=3.8 numpy setuptools scipy tensorflow networkx=2.6 scikit-learn pandas
ENV PATH /opt/conda/envs/dstg/bin:/opt/conda/bin:$PATH
SHELL ["/bin/bash", "-c"]
RUN source activate dstg && pip install git+https://github.com/Su-informatics-lab/DSTG.git@7374fc7aef1471b1c21c4d903a8821ad20bcd2ca
RUN echo "source activate dstg" >> ~/.bashrc
RUN git clone https://github.com/csangara/DSTG.git && \
    # Since we installed tensorflow 2.x, need to import it as v1 compatible
    grep -rl 'import tensorflow' DSTG/DSTG/* | xargs -i@ sed -i 's/import tensorflow/import tensorflow.compat.v1/g' @ && \
    # Due to version conflicts, need to change Float64 to np.float64
    grep -Irl 'Float64' DSTG/DSTG/utils.py | xargs -i@ sed -i "s/'Float64'/np.float64/g" @ && \
    # Add another line to disable eager execution (something that tensorflow 2.x does)
    sed -i '18i tf.compat.v1.disable_eager_execution()' DSTG/DSTG/train.py && \
    # Subset data using Matrix library instead of standard one (for sparse matrices)
    sed -i '2i library(Matrix)' DSTG/DSTG/convert_data.R && \
    # Some changes to the code necessary for it to run
    # Rename columns so they have the same name (for some reason some genes in lab_data2 starts with "X")
    sed -i "20i\    lab_data2 = lab_data2.set_axis(data_train1.columns.values.tolist(), axis=1)" DSTG/DSTG/utils.py && \
    # In our test data we only have 100 spots, which is less than number of neighbors, so add a test for that
    sed -i  "s/combine)/combine,/" DSTG/DSTG/graph.py &&\
    # Changes so DSTG works with smaller amount of spots
    sed -i "82i\    n_spots = count_list2[1].shape[1]" DSTG/DSTG/graph.py && \
    sed -i "s/k_filter=100/k_filter=100 if n_spots >= 100 else n_spots/" DSTG/DSTG/graph.py && \
    sed -i "s/k=30/k=30 if len(spots2) >= 30 else len(spots2)- 1/" DSTG/DSTG/gutils.py && \
    sed -i "s/num_cc=30/num_cc=30 if counts2.shape[1] >=30 else counts2.shape[1]/" DSTG/DSTG/gutils.py && \
    sed -i "s/dims=range(30)/dims=range(30) if len(spots2) >=30 else range(len(spots2))/" DSTG/DSTG/gutils.py && \
    sed -i "56i\                       k_filter=200 if n_spots >= 200 else n_spots)" DSTG/DSTG/graph.py && \
    sed -i "51i\    n_spots = count_list[1].shape[1]" DSTG/DSTG/graph.py && \
    # Filter out count matrix for genes that have 0 counts across all spots
    sed -i "42i\    count_names = count_names.drop([x for x in count_names.values if x not in features])" DSTG/DSTG/gutils.py && \
    # DSTG internally makes use of one function from SPOTlight, which we copy below
    sed -i "1i source('test_spot_fun.R')" DSTG/DSTG/R_utils.R && \
    sed -i 's/SPOTlight:://' DSTG/DSTG/R_utils.R && \
    # Fixing bug for when number of genes <2000
    sed -i "s/return (egen)/return (egen[!is.na(egen)])/" DSTG/DSTG/R_utils.R

# Copy only needed file from SPOTlight instead of installing the whole package
COPY test_spot_fun.R DSTG/DSTG/

WORKDIR "/DSTG/DSTG/"