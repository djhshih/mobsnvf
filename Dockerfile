FROM ubuntu

WORKDIR /home/ubuntu

# Install OS packages
RUN apt-get update && apt-get install -y \
    wget bzip2 \
    r-base r-base-dev \
    gcc git \
    liblzma-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages via CRAN
RUN R -e 'install.packages(c("argparser","stringr"), repos="https://cloud.r-project.org")'

# Install Miniconda
ENV CONDA_DIR=/opt/conda

RUN wget -qO /tmp/mc.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash /tmp/mc.sh -b -p $CONDA_DIR && \
    rm /tmp/mc.sh

ENV PATH="$CONDA_DIR/bin:$PATH"

RUN printf '%s\n' \
    'source /opt/conda/etc/profile.d/conda.sh' \
    'conda activate base' \
    > /etc/profile.d/conda-base.sh

RUN echo -e '\nsource /opt/conda/etc/profile.d/conda.sh && conda activate base' >> /root/.bashrc

# Install GATK4 and gsutil
RUN conda install -y -c bioconda -c conda-forge gatk4 gsutil \
    && conda clean -ya



