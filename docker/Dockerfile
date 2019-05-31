FROM continuumio/miniconda

RUN apt-get -y install ttf-dejavu

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a
