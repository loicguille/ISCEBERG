#Dockerfile

#FROM rocker/shiny:4.1.2


FROM rocker/shiny:4.1.2

RUN apt-get update && apt-get install -y libcurl4-gnutls-dev libssl-dev libxml2 libxml2-dev
RUN apt-get install -y libhdf5-dev
RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'devtools::install_version("ggplot2", version = "3.3.5")'
RUN R -e 'devtools::install_version("Seurat", version = "4.1.0")'
RUN R -e 'devtools::install_version("DT", version = "0.21")'
RUN R -e 'devtools::install_version("waiter", version = "0.2.5")'
RUN R -e 'devtools::install_version("hdf5r", version = "1.3.2")'
RUN R -e 'devtools::install_version("shinycssloaders", version = "1.0.0")'
RUN R -e 'devtools::install_version("clustree", version = "0.4.4")'
RUN R -e 'devtools::install_version("Rcpp", version = "1.0.8")'
RUN R -e 'devtools::install_version("shinymanager", version = "1.0.400")'
RUN R -e 'devtools::install_version("cowplot", version = "1.1.1")'
RUN R -e 'devtools::install_version("RColorBrewer",version = "1.1-2")'
RUN R -e 'devtools::install_version("plotly", version = "4.10.0")'
RUN R -e 'devtools::install_version("dplyr", version = "1.0.8")'
RUN R -e 'devtools::install_version("Matrix", version = "1-4.0")'
RUN R -e 'devtools::install_version("shinyjs", version = "2.1.0")'
RUN R -e 'devtools::install_version("shinyalert", version = "3.0.0")'
RUN R -e 'devtools::install_version("bslib", version = "0.3.1")'
RUN R -e 'devtools::install_version("SeuratObject", version = "4.0.4")'
RUN R -e 'devtools::install_version("colourpicker", version = "1.1.1")'
RUN R -e 'devtools::install_version("rlang",version = "1.0.2")'
RUN R -e 'devtools::install_version("rmarkdown", version = "2.12")'
RUN R -e 'install.packages("igraph", repos = "https://pbil.univ-lyon1.fr/CRAN/")'
RUN R -e 'devtools::install_version("scales", version = "1.1.1")'
RUN R -e 'devtools::install_version("gridExtra", version = "2.3")'
RUN R -e 'devtools::install_version("ggraph", version = "2.0.5")'


COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY /application/ /srv/shiny-server

CMD ["/usr/bin/shiny-server"]
