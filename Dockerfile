FROM rocker/rstudio:3.6.2

RUN mkdir -p /usr/share/man/man1
RUN apt-get update && apt-get install -y libnetcdf-dev libxml2-dev libssl-dev r-cran-rjava=0.9-10-2+b1 default-jdk openbabel=2.4.1+dfsg-3 libopenbabel-dev && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('optparse')"
RUN R -e "BiocManager::install('ChemmineOB')"
RUN R -e "install.packages('string')"
RUN R -e "install.packages('rcdk')"
RUN R -e "install.packages('xml2')"
RUN R -e "install.packages('httr')"
RUN R -e "install.packages('rvest')"
RUN R -e "install.packages('RCurl')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/webchem/webchem_0.2.tar.gz', repos=NULL, type='source')"

WORKDIR /chem_identifier_converter
COPY chem_identifier_converter.r /chem_identifier_converter

