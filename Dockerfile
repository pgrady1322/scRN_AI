########################################
#  Single-Cell Toolkit ─ Dockerfile
#  Build:  docker build -t sc_toolkit:0.1 .
########################################
# ----------  Stage 1  base OS ----------
FROM ubuntu:24.04 AS base
LABEL maintainer="you@example.com"

# Basic build & runtime utilities
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential git curl ca-certificates \
        libgl1  # needed for matplotlib’s Qt back-end in some environments
# Clean up APT cache
RUN rm -rf /var/lib/apt/lists/*

# ----------  Stage 2  Micromamba + env ----------
ARG MAMBA_VER=latest           # 'latest' is a stable alias on the API
ARG MAMBA_ROOT=/opt/conda      # where the env will be stored
ENV PATH=${MAMBA_ROOT}/bin:$PATH
ENV MAMBA_ROOT_PREFIX=${MAMBA_ROOT}

# Install micromamba (download tarball → extract the binary)
RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VER} \
        -o /tmp/micromamba.tar.bz2 && \
    tar -xjf /tmp/micromamba.tar.bz2 -C /usr/local/bin --strip-components=1 bin/micromamba && \
    chmod +x /usr/local/bin/micromamba && \
    rm /tmp/micromamba.tar.bz2

# Copy conda environment spec and create env
COPY env.yml /tmp/env.yml
RUN micromamba create -y -n sc_toolkit --file /tmp/env.yml && \
    micromamba clean -a -y
# Activate by default
ENV CONDA_DEFAULT_ENV=sc_toolkit
ENV PATH=${MAMBA_ROOT}/envs/sc_toolkit/bin:$PATH

# ----------  BLTSA (R) ----------
RUN Rscript -e "options(repos='https://cloud.r-project.org'); install.packages(c('Matrix','FNN','RSpectra','igraph'))" \
    && Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && Rscript -e "BiocManager::install('destiny')" \
    && git clone --depth 1 https://github.com/LiminLi-xjtu/BLTSA.git /opt/BLTSA

# ----------  Stage 3  workflows & CLI ----------
WORKDIR /opt/sc_toolkit

# Copy the source tree
COPY sc_toolkit/ ./sc_toolkit/

# Make it importable without installing
ENV PYTHONPATH=/opt/sc_toolkit:$PYTHONPATH

# RUN pip install --no-cache-dir -e ./sc_toolkit


ENTRYPOINT ["sc_toolkit"]     # main CLI defined below
CMD ["--help"]                # tells users what to do if they just run the image