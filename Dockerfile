# Docker container with conda
FROM continuumio/miniconda3:4.10.3

# Set up python conda environment
COPY .install/python_environment.yml /install/python_environment.yml
RUN conda env create -f /install/python_environment.yml

# Set up R conda environment
COPY .install/R_environment.yml /install/R_environment.yml
RUN conda env create -f /install/R_environment.yml

# Set up R jupyter kernel and make it visible to python
ENV PATH="/opt/conda/envs/py/bin:$PATH"
RUN /opt/conda/envs/R/bin/R -s -e "IRkernel::installspec(sys_prefix = T)"
# Set jupyter data dir for discovering kernels
ENV JUPYTER_DATA_DIR="/opt/conda/envs/py/share/jupyter"

# Make R visible to python environment
ENV PATH="$PATH:/opt/conda/envs/R/bin"

# Make 'py' as default conda environment
RUN sed -i 's/conda activate base/conda activate py/' /root/.bashrc

# Install system dependencies
RUN apt update -y && \ 
    apt install -y git \ 
    curl \ 
    bedtools \ 
    hdf5-tools \ 
    pigz \ 
    pandoc \ 
    graphviz \ 
    zip \ 
    unzip \ 
    libxml2-dev \ 
    libssl-dev \ 
    libcurl4-openssl-dev \ 
    libgit2-dev

# Install hisat2 and bowtie for alignment
RUN conda create -n hisat2 -c bioconda hisat2=2.1.0
RUN conda create -n bowtie -c bioconda bowtie=1.2.2
# to avoid bowtie not finding the right library
# see https://forum.biobakery.org/t/workflow-conda-install-bowtie-2-issue-bowtie2-align-s-error-while-loading-shared-libraries-libtbb-so-2/1831/2
RUN conda install -n bowtie tbb=2020.2

# Make hisat2 and bowtie visible to python environment
ENV PATH="$PATH:/opt/conda/envs/hisat2/bin:/opt/conda/envs/bowtie/bin"