FROM ubuntu

WORKDIR /home/ubuntu

# Install OS packages
RUN apt update -y && apt upgrade -y \
    && apt install -y \
    wget bzip2 git cmake gcc g++ \
    liblzma-dev libcurl4-openssl-dev \
    zlib1g-dev libbz2-dev \
    r-base r-base-dev \
    && rm -rf /var/cache/apt/archives /var/lib/apt/lists/*

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



