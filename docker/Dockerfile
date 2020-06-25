FROM continuumio/miniconda3
SHELL ["/bin/bash", "-c"]


RUN mkdir /pipeline-chromatin-accessibility

ADD ./environment.sh /pipeline-chromatin-accessibility/
RUN bash /pipeline-chromatin-accessibility/environment.sh

ADD ./requirements_pip.txt  /pipeline-chromatin-accessibility/
RUN pip install -r /pipeline-chromatin-accessibility/requirements_pip.txt

ADD ./requirements_R.txt  /pipeline-chromatin-accessibility/
RUN chmod 755 /pipeline-chromatin-accessibility/requirements_R.txt
RUN  /pipeline-chromatin-accessibility/requirements_R.txt


#RUN echo "source activate ATACseq" > ~/.bashrc
#ENV PATH /opt/conda/envs/ATACseq/bin:$PATH
#ENV CONDA_DEFAULT_ENV ATACseq

RUN mkdir /data

WORKDIR /data
