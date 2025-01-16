FROM fedora:40

LABEL org.opencontainers.image.source="https://github.com/marco-compiler/marco-runtime"

ARG LLVM_PARALLEL_COMPILE_JOBS=4
ARG LLVM_PARALLEL_LINK_JOBS=1
ARG LLVM_BUILD_TYPE=Release
ARG LLVM_ENABLE_ASSERTIONS=OFF

# Install compilation dependencies.
RUN dnf update -y && \
    dnf install -y gcc gcc-c++ gfortran perl ninja-build mold cmake ccache \
    git python3-pip python3-virtualenv libxml2-devel wget doxygen

# Create a Python virtual environment.
COPY ./setup_venv.sh /tmp/
RUN chmod +x /tmp/setup_venv.sh && /tmp/setup_venv.sh

# Install LLVM.
COPY ./version_llvm.txt /tmp/
COPY ./install_llvm.sh /tmp/

RUN chmod +x /tmp/install_llvm.sh && \
    cd /root && \
    LLVM_COMMIT=$(cat /tmp/version_llvm.txt) \
    /tmp/install_llvm.sh

# Install solver libraries.
RUN dnf update -y && \
    dnf install -y openblas-devel-0.3.26-4.fc40 \
    suitesparse64-devel-7.6.0-1.fc40 \
    sundials-devel-6.6.2-7.fc40
