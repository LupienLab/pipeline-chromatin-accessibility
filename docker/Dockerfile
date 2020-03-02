FROM continuumio/miniconda3
SHELL ["/bin/bash", "-c"]


RUN mkdir /pipeline-chromatin-accessibility

ADD ./environment.sh /pipeline-chromatin-accessibility/environment.sh
ADD ./initialize.sh /initialize.sh

RUN bash /pipeline-chromatin-accessibility/environment.sh

RUN echo "source activate ATACseq" > ~/.bashrc
ENV PATH /opt/conda/envs/ATACseq/bin:$PATH
ENV CONDA_DEFAULT_ENV ATACseq

RUN mkdir /data


WORKDIR /data
