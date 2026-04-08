FROM rocker/r-ver:4.5.3

RUN /rocker_scripts/setup_R.sh https://packagemanager.posit.co/cran/__linux__/jammy/2026-03-31

# Install system dependencies + BLAST
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install remotes so we can install from GitHub
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')"

# Install your GitHub package
RUN R -e "remotes::install_github('selkamand/blastparse')"
