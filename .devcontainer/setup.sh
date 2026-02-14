#!/bin/bash
# Dev Container setup for Debian packaging with AI support

set -e

echo "🚀 Installing Debian packaging tools..."

apt-get update
apt-get install -y \
    build-essential \
    debhelper \
    devscripts \
    dh-make \
    git-buildpackage \
    pristine-tar \
    pkg-config \
    lintian \
    fakeroot \
    equivs \
    cmake \
    doxygen \
    git

# Install wsclean-specific build dependencies
apt-get install -y \
    libfftw3-dev \
    casacore-dev \
    libgsl-dev \
    libboost-all-dev \
    libcfitsio-dev \
    libhdf5-dev \
    liblapack-dev \
    libopenmpi-dev \
    python3-dev \
    atool

echo "✅ Debian packaging environment initialized!"
echo "Available commands:"
echo "  • debian/fetch-upstream.sh <version>    # Fetch and prepare upstream with submodules"
echo "  • debuild                               # Build the package"
echo "  • lintian <package>.deb                 # Check for packaging issues"
