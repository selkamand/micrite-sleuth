FROM alpine:3.20

RUN apk add --no-cache \
    python3 \
    py3-pip \ 
    bash \
 && python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
 && /opt/venv/bin/pip install --no-cache-dir xopen

ENV PATH="/opt/venv/bin:${PATH}"
ENV PATH="/app:${PATH}"

WORKDIR /app

RUN wget -O index_fasta.py "https://raw.githubusercontent.com/genotoul-bioinfo/dgenies/v1.5.0/src/dgenies/bin/index.py"

RUN chmod +x /app/index_fasta.py

