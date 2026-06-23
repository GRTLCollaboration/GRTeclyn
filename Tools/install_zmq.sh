#!/usr/bin/env bash
# Install libzmq + cppzmq into GRTeclyn/local/zeromq (no root required).
set -euo pipefail

GRTECLYN_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
ZMQ_PREFIX="${GRTECLYN_HOME}/local/zeromq"
BUILD_DIR="${GRTECLYN_HOME}/local/build"
VERSION=4.3.5

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

if [[ ! -f "zeromq-${VERSION}.tar.gz" ]]; then
  wget -q "https://github.com/zeromq/libzmq/releases/download/v${VERSION}/zeromq-${VERSION}.tar.gz"
fi

rm -rf "zeromq-${VERSION}"
tar xzf "zeromq-${VERSION}.tar.gz"
cmake -S "zeromq-${VERSION}" -B "zeromq-${VERSION}/build" \
  -DCMAKE_INSTALL_PREFIX="${ZMQ_PREFIX}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED=ON \
  -DBUILD_STATIC=OFF \
  -DWITH_DOCS=OFF \
  -DENABLE_CURVE=OFF
cmake --build "zeromq-${VERSION}/build" -j"$(nproc)"
cmake --install "zeromq-${VERSION}/build"

wget -q -O "${ZMQ_PREFIX}/include/zmq.hpp" \
  https://raw.githubusercontent.com/zeromq/cppzmq/v4.10.0/zmq.hpp

echo "ZeroMQ installed under ${ZMQ_PREFIX}"
echo "Build RL binary with:"
echo "  cd Examples/RadialRecipe"
echo "  export LD_LIBRARY_PATH=${ZMQ_PREFIX}/lib:\$LD_LIBRARY_PATH"
echo "  make COMP=gnu USE_CUDA=TRUE USE_MPI=FALSE USE_RL=TRUE CUDA_ARCH=90 NVCC_CCBIN=/usr/bin/g++ -j\$(nproc)"
