#!/bin/bash
# Build from clean extracted tarball
# Usage: ./debian/build-from-tarball.sh <VERSION> [DISTRO]

set -e

VERSION="${1}"
DISTRO="${2:-noble}"

if [ -z "$VERSION" ]; then
    echo "Usage: $0 <VERSION> [DISTRO]"
    exit 1
fi

PKG_NAME="wsclean"
TARBALL="../${PKG_NAME}_${VERSION}.orig.tar.gz"
BUILD_DIR="../${PKG_NAME}-${VERSION}"

if [ ! -f "$TARBALL" ]; then
    echo "Error: Tarball not found: $TARBALL"
    exit 1
fi

# Clean previous build directory
if [ -d "$BUILD_DIR" ]; then
    echo "Removing previous build directory..."
    rm -rf "$BUILD_DIR"
fi

# Extract tarball
echo "Extracting tarball..."
cd ..
tar -xzf "${PKG_NAME}_${VERSION}.orig.tar.gz"

# Copy debian directory
echo "Copying debian/ directory..."
cp -r wsclean/debian "${PKG_NAME}-${VERSION}/"

# Build
cd "${PKG_NAME}-${VERSION}"
echo "Building package..."
debuild -us -uc

echo ""
echo "Build complete! Check ../ for .deb files"
