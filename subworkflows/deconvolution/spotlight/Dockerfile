# Build SPOTlight
FROM csangara/seurat:4.1.0
RUN R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-03-17')); remotes::install_github('MarcElosua/SPOTlight', ref = 'spotlight-0.1.7', upgrade='never')"