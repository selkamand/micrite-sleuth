FROM alpine:3.20

RUN apk add --no-cache \
  bash \
  perl \
  curl \
  zlib \
  bzip2 \
  libstdc++ \
  libgcc \
  libgomp \
  ca-certificates \
  coreutils \
  findutils \
  gzip \
  bzip2 \
  tar \
  rsync \
  libc6-compat \
  gcompat \
  build-base \
  zlib-dev \
  bzip2-dev \
  xz-dev \
  perl \
  python3 \
  py3-biopython \
  wget \
  unzip

# Install krakentools
ENV KRAKENTOOLS_VERSION=1.2.1
RUN wget -O "krakentools.zip" "https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v${KRAKENTOOLS_VERSION}.zip" \
  && unzip "krakentools.zip" \
  && cd "/KrakenTools-${KRAKENTOOLS_VERSION}/" \
  && cp *.py /usr/local/bin/


# Install samtools
ENV SAMTOOLS_VERSION=1.22.1
ENV HTSLIB_VERSION=1.22.1

RUN apk add build-base build-base zlib-dev bzip2-dev xz-dev ncurses-dev libcurl && \
  wget -O samtools-${SAMTOOLS_VERSION}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION}/ \
  && ./configure --prefix=/usr/local \
  && make \
  && make install


RUN wget -O htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
  && tar jxvf htslib-${HTSLIB_VERSION}.tar.bz2 \
  && cd htslib-${HTSLIB_VERSION} \
  && ./configure --prefix=/usr/local \
  && make \
  && make install 


# Install seqkit
ENV SEQKIT_VERSION="2.11.0"
RUN wget "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz" \
  && tar xzf "seqkit_linux_amd64.tar.gz" \
  && mv seqkit /usr/local/bin/

# Install Spades & add to path
ENV SPADES_VERSION="4.2.0"
RUN wget "https://github.com/ablab/spades/releases/download/v${SPADES_VERSION}/SPAdes-4.2.0-Linux.tar.gz" \
  && tar xzf "SPAdes-${SPADES_VERSION}-Linux.tar.gz"

ENV PATH=/SPAdes-${SPADES_VERSION}-Linux/bin/:$PATH

# Install R 
RUN apk add 'R=4.4.0-r0'

# Install gnu parallel
RUN apk add parallel

# Install fastqc
ENV FASTQC_VERSION=0.12.1
RUN apk add --no-cache openjdk21 fontconfig ttf-dejavu
RUN wget "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" \
  && unzip "fastqc_v${FASTQC_VERSION}.zip" \
  && cp -r /FastQC/* /usr/local/bin

# Install bwa-mem2
ENV BWAMEM2_VERSION=2.3
RUN wget "https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.3/bwa-mem2-2.3_x64-linux.tar.bz2" \
  && tar xf "bwa-mem2-${BWAMEM2_VERSION}_x64-linux.tar.bz2" \
  && cp /bwa-mem2-${BWAMEM2_VERSION}_x64-linux/* /usr/local/bin/

# Install minimap2
ENV MINIMAP2_VERSION=2.30
RUN wget "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2" \
  &&  tar -jxvf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 \
  && cp minimap2-${MINIMAP2_VERSION}_x64-linux/* /usr/local/bin/


# Install fq
ENV FQ_VERSION=0.12.0
RUN wget "https://github.com/stjude-rust-labs/fq/releases/download/v${FQ_VERSION}/fq-${FQ_VERSION}-x86_64-unknown-linux-gnu.tar.gz" \
  && tar xzf "fq-${FQ_VERSION}-x86_64-unknown-linux-gnu.tar.gz" \
  && cp fq-${FQ_VERSION}-x86_64-unknown-linux-gnu/fq /usr/local/bin

# Install mosdepth
ENV MOSDEPTH_VERSION=0.3.12
RUN wget -O /usr/local/bin/mosdepth "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth" \
  && chmod +x /usr/local/bin/mosdepth

# Install picard
ENV PICARD_VERSION=3.4.0
RUN wget "https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar" \
  && echo echo '#!/usr/bin/env bash' > /usr/local/bin/picard && echo 'java -jar /picard.jar "$@"' >> /usr/local/bin/picard \
  && chmod +x /usr/local/bin/picard

