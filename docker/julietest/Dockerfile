FROM continuumio/miniconda3

LABEL  Eichler Lab

COPY ./environment.yml /

RUN conda env create -f /environment.yml && conda clean -a

ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/conda/envs/assembly_eval/lib/

ENV PATH=/opt/conda/envs/assembly_eval/bin/:${PATH}

RUN apt-get -y update && apt-get -y install file && apt-get install -y gawk && ln -sf /usr/bin/gawk /usr/bin/awk

RUN echo "source activate assembly_eval" > ~/.bashrc

RUN mkdir -p /app/

COPY ./software.list /app/