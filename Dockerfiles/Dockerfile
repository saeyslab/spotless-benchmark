FROM csangara/synthvisium:pre
RUN apt-get install -y git && eval "$(ssh-agent -s)"
RUN R -e "install.packages(c('remotes', 'BiocManager', 'R.utils', 'credentials')); BiocManager::install('DropletUtils')"
# Run it once inside the container and retrieve the known_hosts file
COPY known_hosts /root/.ssh/ 
RUN R -e "library(credentials); remotes::install_git('git@github.com:csangara/synthvisium.git', upgrade='never')"

# First step - create csangara/synthvisium:pre
# FROM csangara/seurat:3.2.3
# RUN ssh-keygen -t rsa -f ~/.ssh/id_rsa
# Then add the public key to deploy key in csangara/synthvisium