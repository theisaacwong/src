FROM ubuntu:22.04

WORKDIR /home

#### Basic image utilities
RUN apt-get update && \
    apt-get full-upgrade -y && \
    apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    wget \
    curl \
    bc \
    unzip \
    bzip2 \
    less \
    bedtools \
    bcftools \
    gcc \
    tabix \
    jq \
    git \
    gpg-agent \
    build-essential \
    openjdk-17-jdk \
    cmake \
    clang \
    libclang-dev \
    libncurses-dev \
    zlib1g-dev \
    libbz2-dev \
    libssl-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    make \
    cmake \
    vim \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

# Install R
ARG R_VERSION=4.3.3
ARG OS_IDENTIFIER=ubuntu-2204

RUN wget https://cdn.posit.co/r/${OS_IDENTIFIER}/pkgs/r-${R_VERSION}_1_amd64.deb && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -f -y ./r-${R_VERSION}_1_amd64.deb && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/${R_VERSION}/lib/R /usr/lib/R && \
    rm r-${R_VERSION}_1_amd64.deb && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('parallel',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install miniconda
ENV DOWNLOAD_DIR /downloads
ENV CONDA_URL https://repo.anaconda.com/miniconda/Miniconda3-py310_23.10.0-1-Linux-x86_64.sh
ENV CONDA_SHA256 "c7a34df472feb69805b64df6e8db58363c5ccab41cd3b40b07e3e6dfb924359a"
ENV CONDA_PATH /opt/miniconda
ENV PATH $CONDA_PATH/bin:$PATH
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "$(sha256sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS=' ' -v ORS='' '{print $1}')" = "$CONDA_SHA256" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh && \
    conda clean -afy && \
    conda config --set auto_update_conda false && \
    conda config --set solver libmamba && \
    rm -rf /root/.cache/pip

RUN wget  https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -xvf bwa-0.7.17.tar.bz2 && \
    rm -rf bwa-0.7.17.tar.bz2 && \
    ln -s /home/bwa-0.7.17/bwakit/run-bwamem /home/bwa-0.7.17/bwakit/bwa
ENV PATH ./bwa-0.7.17/bwakit/:$PATH

# Install htslib
WORKDIR /opt/htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 && \
    tar -xvf htslib-1.19.1.tar.bz2 && \   
    rm -rf htslib-1.19.1.tar.bz2
WORKDIR /opt/htslib/htslib-1.19.1 
RUN ./configure && \
    make 

# Install samtools
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    tar -xvf samtools-1.19.2.tar.bz2 && \
    rm -rf samtools-1.19.2.tar.bz2
WORKDIR /opt/samtools/samtools-1.19.2
RUN ./configure --with-htslib=/opt/htslib/htslib-1.19.1 && \
    make && \
    make install && \
    mkdir -p /home/cache
ENV REF_PATH="/opt/samtools/samtools-1.19.2/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
ENV REF_CACHE="/opt/samtools/samtools-1.19.2/ref/cache/%2s/%2s/%s"

COPY Dockerfile .
