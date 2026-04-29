FROM continuumio/miniconda3@sha256:42cd2ca0ece04579b6127e1a1a0f83c25a079145d408eb65e39206ac05069a77
LABEL tobias 0.16.0 uropa 4.0.2

ENV SHELL=bash

RUN mkdir -p /opt/lib
WORKDIR /opt/lib


RUN apt-get update &&  apt-get install -y --no-install-recommends autoconf automake bzip2 dpkg-dev file g++ gcc imagemagick libbz2-dev libc6-dev libcurl4-openssl-dev libdb-dev libevent-dev libffi-dev libgdbm-dev libgeoip-dev libglib2.0-dev libgmp-dev libjpeg-dev libkrb5-dev liblzma-dev libmagickcore-dev libmagickwand-dev libncurses5-dev libncursesw5-dev libpng-dev libpq-dev libreadline-dev libsqlite3-dev libssl-dev libtool libwebp-dev libxml2-dev libxslt-dev libyaml-dev make patch unzip xz-utils zlib1g-dev  && rm -rf /var/lib/apt/lists/*

ENV PATH /opt/conda/envs/env/bin:$PATH
COPY environment.yml /opt/lib
RUN conda env create -f /opt/lib/environment.yml



RUN echo "source activate env" > ~/.bashrc

WORKDIR /data

