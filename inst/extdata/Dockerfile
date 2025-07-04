FROM rocker/shiny:4.3.2

RUN apt-get update \
        && apt-get install -y \
        apt-utils \
        manpages-dev \
        libnetcdf-dev \
        libxml2-dev \
        libglpk-dev \
        librsvg2-dev \
        libopenbabel-dev \
        libnode-dev  \
        libz-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libxt-dev \
        libglpk-dev \
        libgmp-dev \
        libboost-all-dev \
        g++ \
        && apt-get clean


# Copy the script into the container
COPY babel_installation.sh /usr/local/bin/babel_installation.sh

# Make it executable
RUN chmod +x /usr/local/bin/babel_installation.sh

RUN R -e "install.packages('BiocManager')"
# RUN R -e "BiocManager::install('omicsViewer')"
RUN R -e "BiocManager::install('xcms')"
RUN R -e "BiocManager::install('MAIT')"
RUN R -e "install.packages('devtools')"
RUN R -e 'install.packages("remotes")'

RUN R -e 'remotes::install_version("survminer", version = "0.4.9")'
RUN R -e 'remotes::install_version("survival", version = "3.3-1")'
RUN R -e 'remotes::install_version("fastmatch", version = "1.1-3")'
RUN R -e 'remotes::install_version("reshape2", version = "1.4.4")'
RUN R -e 'remotes::install_version("beeswarm", version = "0.4.0")'
# grDevices is part of base R — do NOT install it
RUN R -e 'remotes::install_version("shinycssloaders", version = "1.0.0")'
RUN R -e 'remotes::install_version("shinythemes", version = "1.2.0")'
RUN R -e 'remotes::install_version("networkD3", version = "0.4")'
RUN R -e 'remotes::install_version("httr", version = "1.4.3")'
RUN R -e 'remotes::install_version("RColorBrewer", version = "1.1-3")'
RUN R -e 'remotes::install_version("psych", version = "2.2.5")'
RUN R -e 'remotes::install_version("stringr", version = "1.4.0")'
RUN R -e 'remotes::install_version("shiny", version = "1.7.1")'
RUN R -e 'remotes::install_version("shinydashboard", version = "0.7.2")'
RUN R -e 'remotes::install_version("shinyWidgets", version = "0.7.0")'
RUN R -e 'remotes::install_version("shinybusy", version = "0.3.1")'
RUN R -e 'remotes::install_version("matrixStats", version = "0.62.0")'
RUN R -e 'remotes::install_version("flatxml", version = "0.1.1")'
RUN R -e 'remotes::install_version("excelR", version = "0.4.0")'
RUN R -e 'remotes::install_version("shinyjs", version = "2.1.0")'
RUN R -e 'remotes::install_version("shinyFiles", version = "0.9.2")'
RUN R -e 'remotes::install_version("DT", version = "0.23")'
RUN R -e 'remotes::install_version("plotly", version = "4.10.0")'
RUN R -e 'remotes::install_version("openxlsx", version = "4.2.5")'
RUN R -e 'remotes::install_version("yaml", version = "2.3.5")'
RUN R -e 'remotes::install_version("curl", version = "4.3.2")'
RUN R -e 'remotes::install_version("sortable", version = "0.4.5")'
RUN R -e 'remotes::install_version("BiocManager", version = "1.30.18")'
RUN R -e 'remotes::install_version("password", version = "1.0-0")'
RUN R -e 'remotes::install_version("ggseqlogo", version = "0.1")'
RUN R -e 'remotes::install_version("devtools", version = "2.4.3")'
RUN R -e 'remotes::install_version("Matrix", version = "1.4-1")'
RUN R -e 'remotes::install_version("RSQLite", version = "2.2.14")'

# Install BiocManager and fgsea (from Bioconductor 3.17)
RUN R -e  "BiocManager::install('fgsea')"

RUN /usr/local/bin/babel_installation.sh
RUN R -e 'BiocManager::install(c("Biobase",  "S4Vectors", "SummarizedExperiment","ChemmineR","ChemmineOB"), update = FALSE)'
RUN R -e 'devtools::install_github("mengchen18/omicsViewer", dependencies = TRUE)'
RUN R -e "devtools::install_github('nava20ir/xcmsViewer', upgrade='never')"



RUN mkdir /home/shiny/data && \
        mkdir /home/shiny/app

#COPY *.R /home/shiny/app/
# Copy your app
COPY *.R /srv/shiny-server/
RUN sed -i 's/3838/5968/' /etc/shiny-server/shiny-server.conf

EXPOSE 5968

# Run shiny server
CMD ["/usr/bin/shiny-server"]
