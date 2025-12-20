FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    autotools-dev \
    libtool \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    bzip2 \
    make \
    sudo \
    git \
    wget \
    ca-certificates \
    libboost-program-options-dev \
    libboost-iostreams-dev \
    libboost-serialization-dev \
    && rm -rf /var/lib/apt/lists/*

# Build HTSlib
WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -xjf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    ./configure --enable-libcurl && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd .. && rm -rf htslib-1.16*


# Install cmake
RUN mkdir -p /usr/local
RUN cd /usr/local &&\
    wget https://github.com/Kitware/CMake/releases/download/v4.0.2/cmake-4.0.2-linux-x86_64.sh &&\
    chmod +x cmake-4.0.2-linux-x86_64.sh &&\
    ./cmake-4.0.2-linux-x86_64.sh --skip-license

# Install Doxygen 1.9.5
RUN wget https://www.doxygen.nl/files/doxygen-1.9.5.linux.bin.tar.gz && \
    tar -xzf doxygen-1.9.5.linux.bin.tar.gz && \
    ln -sf /opt/doxygen-1.9.5/bin/doxygen /usr/local/bin/doxygen && \
    rm doxygen-1.9.5.linux.bin.tar.gz


# Create vscode user
RUN useradd -m -s /bin/bash vscode && \
    echo "vscode ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/vscode && \
    chmod 0440 /etc/sudoers.d/vscode

USER vscode
WORKDIR /workspace