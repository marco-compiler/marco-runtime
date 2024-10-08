FROM debian:12

LABEL org.opencontainers.image.source="https://github.com/marco-compiler/marco-runtime"

ARG LLVM_PARALLEL_COMPILE_JOBS=4
ARG LLVM_PARALLEL_LINK_JOBS=1

RUN apt update -y && \
    apt install -y build-essential gfortran ninja-build lld mold cmake ccache \
    git python3-pip python3-venv libxml2-dev libtinfo-dev wget doxygen

COPY ./setup_venv.sh /tmp/
RUN chmod +x /tmp/setup_venv.sh && /tmp/setup_venv.sh

COPY ./version_llvm.txt /tmp/
COPY ./install_llvm.sh /tmp/

RUN chmod +x /tmp/install_llvm.sh && \
    cd /root && \
    LLVM_COMMIT=$(cat /tmp/version_llvm.txt) \
    LLVM_BUILD_TYPE=Release \
    LLVM_ENABLE_ASSERTIONS=ON \
    /tmp/install_llvm.sh

RUN apt update -y && \
    apt install -y libopenblas-dev=0.3.21+ds-4 \
    libsuitesparse-dev=1:5.12.0+dfsg-2 \
    libsundials-dev=6.4.1+dfsg1-3
