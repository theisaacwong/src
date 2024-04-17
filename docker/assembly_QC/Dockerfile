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
    samtools \
    bcftools \
    gcc \
    tabix \
    jq \
    git \
    gpg-agent \
    build-essential \
    openjdk-17-jdk \
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
RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('parallel',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('scales',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('argparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

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

# Install R packages via conda, these are mainly bioconductor packages
RUN conda install -c bioconda -c r -c conda-forge bioconductor-genomicranges r-gtools r-ggplot2 r-ggnewscale

# Install blast
RUN wget -nv https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN tar zxvpf ncbi-blast-2.15.0+-x64-linux.tar.gz
ENV PATH="${PATH}:/home/ncbi-blast-2.15.0+/bin"

#  Install Meryl and winnowmap
RUN conda install -c bioconda -y winnowmap=2.03 \
    minimap2=2.28-0 \
    meryl=1.4.1
#RUN wget https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz
#RUN tar -xJf meryl-1.4.1.Linux-amd64.tar.xz
#ENV PATH="${PATH}:/home/meryl-1.4.1/bin"

# Install pcre2 and hdf5. note: might need to install pip versions instead if these don't work
RUN apt-get update && \ 
    apt-get install -y pcre2-utils libpcre2-dev hdf5-tools

# Install IGV
RUN wget https://data.broadinstitute.org/igv/projects/downloads/2.17/IGV_2.17.4.zip
RUN unzip IGV_2.17.4.zip
ENV PATH="${PATH}:/home/IGV_2.17.4"

# Install Merqury
RUN wget https://github.com/marbl/merqury/archive/v1.3.tar.gz
RUN tar -zxvf v1.3.tar.gz
ENV PATH="${PATH}:/home/merqury-1.3"


# Install compleasm
RUN conda create -n compleasm -c conda-forge -c bioconda compleasm --yes
RUN conda init bash
# RUN echo "conda activate compleasm" >> ~/.bashrc
ENV PATH=/opt/conda/envs/compleasm/bin:$PATH

# Install rustybam
RUN wget -nv https://github.com/mrvollger/rustybam/releases/download/v0.1.33/rustybam-x86_64-unknown-linux-gnu.tar.xz
RUN tar -xJf rustybam-x86_64-unknown-linux-gnu.tar.xz
ENV PATH=/home/rustybam-x86_64-unknown-linux-gnu:$PATH

RUN rm -rf IGV_2.17.4.zip \
    meryl-1.4.1.Linux-amd64.tar.xz \
    ncbi-blast-2.15.0+-x64-linux.tar.gz \
    v1.3.tar.gz

# miniconda is 4.3.1, want ot change so 4.3.3 is used by default and also can use packages installed by conda
ENV PATH=/opt/R/4.3.3/bin:$PATH
RUN echo ".libPaths(c('/opt/R/4.3.3/lib/R/library', '/opt/miniconda/lib/R/library',.libPaths()))" > /home/.Rprofile

COPY Dockerfile .
