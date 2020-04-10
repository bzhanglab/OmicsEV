FROM bioconductor/bioconductor_docker:RELEASE_3_10


MAINTAINER wenbostar@gmail.com


ARG PYTHON=python3
ARG PIP=pip3

# See http://bugs.python.org/issue19846
ENV LANG C.UTF-8

RUN  rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    unzip \
    pkg-config \
    curl \
    ${PYTHON} \
    ${PYTHON}-pip && \
    apt-get clean

RUN ${PIP} install --upgrade pip
RUN ${PIP} install --upgrade setuptools

RUN ln -s $(which ${PYTHON}) /usr/local/bin/python

RUN ${PIP} install numpy pandas matplotlib seaborn scikit-learn

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get -y  install --fix-missing tcl8.6-dev \
    tk \
    expat \
    libexpat-dev && \
    apt-get clean





ADD install.R /tmp/

RUN R -f /tmp/install.R

RUN echo "R_LIBS=/usr/local/lib/R/host-site-library:\${R_LIBS}" > /usr/local/lib/R/etc/Renviron.site
RUN echo "R_LIBS_USER=''" >> /usr/local/lib/R/etc/Renviron.site
RUN echo "options(defaultPackages=c(getOption('defaultPackages'),'BiocManager'))" >> /usr/local/lib/R/etc/Rprofile.site

RUN chmod -R 755 /opt/

#change working directory
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
