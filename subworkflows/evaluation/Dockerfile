# Build MuSiC
FROM rocker/tidyverse:3.6.3
RUN R -e "install.packages(c('R.utils', 'precrec'))"
ENTRYPOINT [ "R" ]