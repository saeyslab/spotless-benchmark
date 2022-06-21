# Build spatialDWLS
FROM csangara/seurat:4.1.0
COPY requirements.txt /tmp/
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    python3.6 python3-pip python3-setuptools python3-dev libmagick++-dev \
    && pip3 install -r /tmp/requirements.txt
RUN R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-04-21')); install.packages(c('quadprog', 'magick'))"
RUN R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-04-21')); remotes::install_github('RubD/Giotto@eb729dba433abbf9059c442ddff2032769260980', upgrade='never')"