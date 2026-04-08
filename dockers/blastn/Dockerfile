FROM alpine:3.23 AS build

ARG BLAST_VERSION=2.17.0

RUN apk add --no-cache \
    ca-certificates \
    curl \
    perl \
    build-base \
    sqlite-dev \
    lmdb-dev \
    zstd-dev \
    bzip2-dev \
    zlib-dev

WORKDIR /tmp

RUN curl -fsSLO "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${BLAST_VERSION}+-src.tar.gz" \
 && tar -xzf "ncbi-blast-${BLAST_VERSION}+-src.tar.gz" \
 && srcdir="$(find . -maxdepth 1 -type d -name 'ncbi-blast-*' | head -n 1)" \
 && mv "$srcdir" /opt/ncbi-blast \
 && rm -f "ncbi-blast-${BLAST_VERSION}+-src.tar.gz"

WORKDIR /opt/ncbi-blast/c++

RUN ./configure \
 && cd ReleaseMT/build \
 && make all_r

FROM alpine:3.23

RUN apk add --no-cache \
    ca-certificates \
    libstdc++ \
    libgomp \
    sqlite-libs \
    lmdb \
    zstd-libs \
    libbz2 \
    zlib

COPY --from=build /opt/ncbi-blast /opt/ncbi-blast

ENV PATH="/opt/ncbi-blast/c++/ReleaseMT/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/ncbi-blast/c++/ReleaseMT/lib" \
    BLASTDB="/blast/db"

WORKDIR /work

RUN mkdir -p /blast/db /work \
 && blastn -version

ENTRYPOINT ["blastn"]
CMD ["-help"]

