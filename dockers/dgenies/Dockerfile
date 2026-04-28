FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="${VIRTUAL_ENV}/bin:/app:${PATH}"

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-venv \
    python3-pip \
    wget \
    bash \
    build-essential \
    libffi-dev \
    libssl-dev \
    zlib1g-dev \
    libjpeg-dev \
    libfreetype6-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv "${VIRTUAL_ENV}" \
 && "${VIRTUAL_ENV}/bin/pip" install --no-cache-dir --upgrade pip setuptools wheel \
 && "${VIRTUAL_ENV}/bin/pip" install --no-cache-dir xopen \
 && "${VIRTUAL_ENV}/bin/pip" install --no-cache-dir dgenies

WORKDIR /app

# Install dgenies index scripts (used on query and fasta sequences)
RUN wget -O /app/index_fasta.py \
    "https://raw.githubusercontent.com/genotoul-bioinfo/dgenies/v1.5.0/src/dgenies/bin/index.py" \
 && chmod +x /app/index_fasta.py
