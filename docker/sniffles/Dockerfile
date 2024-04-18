FROM ubuntu:22.04

WORKDIR /home

#### Basic image utilities
RUN apt-get update && \
    apt-get full-upgrade -y && \
    apt-get install -y --no-install-recommends \
    wget \
    curl \
    bc \
    unzip \
    bzip2 \
    less \
    g++ \
    tabix \
    jq \
    git \
    build-essential \
    vim \
    cmake \
    clang \
    python3 \
    python3-pip \
    libclang-dev \
    libncurses-dev \
    zlib1g-dev \
    libbz2-dev \
    libssl-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    make \
    cmake \
    autoconf \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

RUN pip install sniffles==2.3.2 pysam psutil



WORKDIR /home
COPY Dockerfile .
