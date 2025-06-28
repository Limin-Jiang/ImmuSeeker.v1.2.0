# Use an official Python runtime as a parent image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update

# Install system dependencies required for systemfonts and ggplot-related packages
RUN apt-get update && apt-get install -y \
	libcairo2-dev  \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgl1-mesa-dev \
    libxt-dev \
    build-essential \
    && apt-get clean

FROM rocker/r-base:4.3.2


RUN R -e "install.packages(c('data.table', 'dplyr'))"
RUN R -e "install.packages(c('docopt'))"
RUN R -e "install.packages(c('stringr'))"
RUN apt-get update && apt-get install -y \
	libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    zlib1g-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev
RUN R -e "install.packages('RCurl', repos = 'https://cloud.r-project.org')"
RUN R -e "install.packages('BiocManager'); BiocManager::install('Biostrings')"
RUN R -e "install.packages(c('base64enc', 'jsonlite', 'cpp11'), repos = 'https://cloud.r-project.org')"
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libpango1.0-dev \
    build-essential
	
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/systemfonts_1.2.3.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages(c('ggforce'), repos='https://cloud.r-project.org')"
RUN R -e "install.packages(c('ggraph'), repos='https://cloud.r-project.org')"
RUN R -e "install.packages('BiocManager'); BiocManager::install('phyloseq')"
RUN R -e "install.packages(c('igraph'))"

RUN apt -y install samtools
RUN apt -y install bowtie

# Set the working directory
WORKDIR /app

RUN mkdir /ImmuSeeker_data
#COPY ./data /ImmuSeeker_data
# Copy the current directory contents into the container at /app
COPY . /app
#COPY dockerfile /app
#COPY ImmuSeeker /app

# Define default command to run when the container starts
#CMD ["bash"]
COPY ImmuSeeker /usr/local/bin/ImmuSeeker
RUN chmod +x /usr/local/bin/ImmuSeeker
RUN chown -R root:root /app && chmod -R 777 /app
#USER root
RUN chown -R root:root /ImmuSeeker_data && chmod -R 777 /ImmuSeeker_data
#CMD ["R"]

ENTRYPOINT ["ImmuSeeker"]
