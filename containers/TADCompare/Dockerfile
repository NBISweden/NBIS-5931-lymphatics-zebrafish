FROM rocker/tidyverse@sha256:795f7ea93c8028a865e25456575db4e02377e35d5e285aeae9b2c0ed4f7ba4a8
LABEL tadcompare_433 based on rocker/tidyverse:4.4.3

WORKDIR .

RUN mkdir ../opt/lib

RUN apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y \
 git \
 git-lfs \
 make \
 gzip \
 rename \
 libglpk-dev \
 libfftw3-dev \
 libxt6 \
 gcc \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN /rocker_scripts/install_tidyverse.sh

#packages from Bioconductor
COPY install_Bioconductorlibs_tadcompare.R ../opt/lib
RUN Rscript ../opt/lib/install_Bioconductorlibs_tadcompare.R

ENV SHELL=bash

WORKDIR /data


