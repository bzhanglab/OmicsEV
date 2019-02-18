FROM bioconductor/release_mscore2


MAINTAINER wenbostar@gmail.com

RUN  rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*


RUN apt-get update && \
    apt-get -y  install --fix-missing tcl8.6-dev \
    tk




ADD install.R /tmp/

RUN R -f /tmp/install.R

RUN echo "R_LIBS=/usr/local/lib/R/host-site-library:\${R_LIBS}" > /usr/local/lib/R/etc/Renviron.site
RUN echo "R_LIBS_USER=''" >> /usr/local/lib/R/etc/Renviron.site
RUN echo "options(defaultPackages=c(getOption('defaultPackages'),'BiocManager'))" >> /usr/local/lib/R/etc/Rprofile.site

