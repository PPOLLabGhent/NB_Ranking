# syntax=docker/dockerfile:1

FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package manager and install required tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    r-base \
    wget \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

FROM python:3.12-slim
WORKDIR /code

RUN apt-get update && apt-get install -y git r-base wget
RUN git clone https://github.com/PPOLLabGhent/lostdata.git
ENV SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
RUN git clone https://github.com/PPOLLabGhent/bidali.git
RUN cd lostdata && pip install --no-cache-dir -e .
RUN cd bidali && pip install --no-cache-dir --no-cache-dir -e .
RUN pip install leopard openpyxl scikit-learn

COPY . /code/nb_ranking

# RUN for R libraries
RUN Rscript -e "install.packages('BiocManager'); \
    BiocManager::install('limma'); \
    BiocManager::install('preprocessCore')"

#rmdir /root/.local/share/lostdata/private
#ln -s /root/LSData /root/.local/share/lostdata/private
