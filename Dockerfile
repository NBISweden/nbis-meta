FROM centos:7.6.1810

LABEL maintainer="John Sundh" email=john.sundh@nbis.se
LABEL description="Dockerfile to reproduce SFB bin from human data"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /analysis

# Update packages
RUN yum -y update && yum -y install bzip2 git && yum clean all

# Install miniconda
RUN curl https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O && \
    bash Miniconda3-4.6.14-Linux-x86_64.sh -bf -p /opt/miniconda3/ && \
    rm Miniconda3-4.6.14-Linux-x86_64.sh && \
    ln -s /opt/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/miniconda3/etc/profile.d/conda.sh " >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Add conda to PATH and set locale
ENV PATH="/opt/miniconda3/bin:${PATH}"
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8

# Add repository files
COPY source ./nbis-meta
COPY environment.yml ./nbis-meta
COPY config ./nbis-meta
COPY config_*.yml ./nbis-meta
COPY gut_metagenome_samples.tsv ./nbis-meta
COPY infant_gut_metagenome_samples.tsv ./nbis-meta
COPY sfb_metagenome_novaseq_preprocessed_samples.tsv ./nbis-meta
COPY sfb_metagenome_novaseq_samples.tsv ./nbis-meta
COPY sfb_metagenome_samples.tsv ./nbis-meta
COPY zhejiang_infant_gut_metagenome_samples.tsv ./nbis-meta
COPY run.sh ./nbis-meta
COPY Snakefile ./nbis-meta

# Install environment
RUN conda env create -f /code/nbis-meta/environment.yml -n nbis-meta && conda clean --all -y
RUN source activate nbis-meta

CMD /bin/bash run.sh