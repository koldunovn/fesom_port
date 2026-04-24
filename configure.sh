#!/usr/bin/env bash
# Login-shell build script. Run with `bash -l configure.sh [--clean]`.
# Mirrors fesom2/configure.sh layout.

set -e

HERE=$PWD
SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
BUILD_DIR=${BUILD_DIR:-build}

CLEAN_BUILD=false
CMAKE_EXTRA=()
for arg in "$@"; do
    if [ "$arg" = "--clean" ]; then
        CLEAN_BUILD=true
    else
        CMAKE_EXTRA+=("$arg")
    fi
done

if [ "$CLEAN_BUILD" = true ]; then
    echo "Cleaning ${BUILD_DIR}..."
    rm -rf "${BUILD_DIR}"
fi

source "${SOURCE_DIR}/env.sh"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

cmake "${SOURCE_DIR}" -DCMAKE_BUILD_TYPE=Release "${CMAKE_EXTRA[@]}"
make -j"$(nproc --all)"

echo
echo "Built: ${SOURCE_DIR}/${BUILD_DIR}/fesom_port"
