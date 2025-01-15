FROM ghcr.io/prefix-dev/pixi:latest

USER root

ARG FIRE_VERSION=v0.1.1

ARG PIXI_VERSION=0.39.5

RUN apt-get update && \
    apt-get full-upgrade -y && \
    apt-get install -y --no-install-recommends \
    wget \
    less \
    curl \
    bc \
    unzip \
    bzip2 \
    gcc \
    jq \
    git \
    gpg-agent \
    build-essential \
    libclang-dev \
    libncurses-dev \
    zlib1g-dev \
    libbz2-dev \
    libssl-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    vim \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/local/bin/pixi:${PATH}"

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
# Install Miniconda on x86 or ARM platforms
RUN arch=$(uname -m) && \
    if [ "$arch" = "x86_64" ]; then \
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"; \
    elif [ "$arch" = "aarch64" ]; then \
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"; \
    else \
    echo "Unsupported architecture: $arch"; \
    exit 1; \
    fi && \
    wget $MINICONDA_URL -O miniconda.sh && \
    mkdir -p /root/.conda && \
    bash miniconda.sh -b -p /root/miniconda3 && \
    rm -f miniconda.sh
RUN conda --version

WORKDIR /opt/conda/envs/snakemake

RUN conda create -c conda-forge -c bioconda -n snakemake 'snakemake>=8.4'

RUN conda env list > /conda_env.txt

# RUN curl -Ls \
#     "https://github.com/prefix-dev/pixi/releases/download/v${PIXI_VERSION}/pixi-$(uname -m)-unknown-linux-musl" \
#     -o /pixi && chmod +x /pixi
# RUN /pixi --version
# ENV PATH="/root/.pixi/bin:${PATH}"

WORKDIR /opt/
RUN git clone https://github.com/fiberseq/FIRE.git
WORKDIR /opt/FIRE 
RUN git checkout ${FIRE_VERSION}

RUN /usr/local/bin/pixi install

RUN mkdir /opt/fire

ENV PATH="/root/miniconda3/envs/snakemake/bin:/root/miniconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/fire"

RUN echo "#!/bin/bash\nmkdir -p \${HOME}/fire_env;\nexport SNAKEMAKE_CONDA_PREFIX=\${HOME}/fire_env;\nexport TMP=\${HOME}/fire_env;\nexport TEMPDIR=\${HOME}/fire_env;\nexport TEMP=\${HOME}/fire_env;\nexport APPTAINER_CACHEDIR=\${HOME}/fire_env;\nexport PATH='/root/miniconda3/envs/snakemake/bin:/root/miniconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/fire:/opt/FIRE';\nconda activate /root/miniconda3/envs/snakemake" > /opt/fire/load_fire.sh

RUN chmod +x /opt/fire/load_fire.sh

ENTRYPOINT [". /opt/fire/load_fire.sh"]

COPY Dockerfile .
