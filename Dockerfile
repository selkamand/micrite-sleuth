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

#\
 # && chmod +x "*" \
 # && cp * /usr/local/bin/


# Add scripts path
# ENV PATH="$PATH:/app"

# WORKDIR /app
