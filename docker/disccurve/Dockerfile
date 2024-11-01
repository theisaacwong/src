FROM continuumio/miniconda3

COPY ./environment.yml /

RUN conda env create -f /environment.yml && conda clean -a

ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/conda/envs/align-basics/lib/

ENV PATH=/opt/conda/envs/align-basics/bin/:${PATH}

RUN echo "source activate disc-curve" > ~/.bashrc

RUN mkdir -p /app/

COPY ./Dockerfile /
