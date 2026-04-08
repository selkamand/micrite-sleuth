# Use a stable R base image (Ubuntu-based)
FROM rocker/r-ver:4.5.3

# Install system dependencies
# 1. ncbi-blast+ for blastn
# 2. libcurl/libssl/libxml2 for R package installation (devtools/remotes)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install 'remotes' package to allow installing from GitHub
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org/')"

# Install your R package from GitHub
RUN R -e "remotes::install_github('selkamand/blastparse')"

# Set the default command to R
CMD ["R"]
