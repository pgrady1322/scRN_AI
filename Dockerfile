########################################
#  Single-Cell Toolkit ─ Dockerfile
#  Build:  docker build -t scrn_ai:0.1 .
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
RUN micromamba create -y -n scrn_ai --file /tmp/env.yml && \
    micromamba clean -a -y
# Activate by default
ENV CONDA_DEFAULT_ENV=scrn_ai
ENV PATH=${MAMBA_ROOT}/envs/scrn_ai/bin:$PATH

# ----------  BLTSA (R) ----------
RUN Rscript -e "options(repos='https://cloud.r-project.org'); install.packages(c('Matrix','FNN','RSpectra','igraph'))" \
    && Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && Rscript -e "BiocManager::install('destiny')" \
    && git clone --depth 1 https://github.com/LiminLi-xjtu/BLTSA.git /opt/BLTSA

# ----------  Stage 3  workflows & CLI ----------
WORKDIR /opt/scrn_ai

# Copy the source tree and setup files
COPY scrn_ai/ ./scrn_ai/
COPY setup.py ./setup.py

# Install the package (creates the scrn_ai command)
RUN pip install --no-cache-dir -e .

ENTRYPOINT ["scrn_ai"]     # main CLI defined below
CMD ["--help"]                # tells users what to do if they just run the image