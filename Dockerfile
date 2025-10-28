# Base image: Miniconda (for Snakemake + bio tools)
FROM continuumio/miniconda3:latest

# Metadata
LABEL maintainer="Dhrupad Banerjee <245hsbb010@ibab.ac.in>"
LABEL description="Dockerized RNA-Seq Snakemake Pipeline"

# Set working directory
WORKDIR /pipeline

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    graphviz \
    && rm -rf /var/lib/apt/lists/*

# Install Snakemake and bioinformatics tools
RUN conda install -y -c bioconda -c conda-forge \
    snakemake \
    fastp \
    hisat2 \
    samtools \
    subread \
    bioconductor-deseq2 \
    r-ggplot2 \
    r-ggrepel \
    r-yaml

# Copy pipeline files into container
COPY . /pipeline

# Set default command
ENTRYPOINT ["snakemake"]
CMD ["--help"]

