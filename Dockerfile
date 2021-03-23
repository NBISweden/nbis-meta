FROM continuumio/miniconda3:4.8.2

LABEL maintainer="John Sundh" email=john.sundh@nbis.se
LABEL description="Docker image for a snakemake workflow for metagenomics"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /analysis

# Set tmpdir
ENV TMPDIR="/scratch"

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && apt-get clean

# Add environment file
COPY environment.yml .

# Install environment into base
RUN conda env update -n base -f environment.yml && conda clean -a

# Add workflow
COPY config config
COPY workflow workflow

# Run workflow
ENTRYPOINT ["snakemake", "--use-conda", "-j", "1"]
