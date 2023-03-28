#Dockerfile

FROM rocker/shiny:4.1.2


RUN apt-get update && apt-get install -y libcurl4-gnutls-dev libssl-dev libxml2
RUN apt-get install -y libhdf5-dev
RUN R -e 'install.packages(c("BiocManager","ggplot2", "Seurat","DT", "waiter","hdf5r", "shinycssloaders", "clustree", "Rcpp", "remotes", "shinymanager", "cowplot", "RColorBrewer","plotly","dplyr","Matrix","shinyjs","shinyalert","bslib","SeuratObject","colourpicker","rlang","rmarkdown"))'
RUN R -e 'remotes::install_github("igraph/rigraph")'

COPY shiny-server.conf /etc/shiny-server/shiny-server.conf	
COPY application/ /srv/shiny-server/

CMD ["/usr/bin/shiny-server"]
