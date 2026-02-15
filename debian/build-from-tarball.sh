#!/bin/bash
# Build from clean extracted tarball
# Usage: ./debian/build-from-tarball.sh <VERSION> [DISTRO]

set -e

VERSION="${1}"
DISTRO="${2:-noble}"
INSTALL_DEPS="${3:-yes}"

if [ -z "$VERSION" ]; then
    echo "Usage: $0 <VERSION> [DISTRO] [INSTALL_DEPS]"
    echo "Example: $0 3.6 noble yes"
    echo ""
    echo "INSTALL_DEPS: yes (default) or no"
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

# Install build dependencies if requested
if [ "$INSTALL_DEPS" = "yes" ]; then
    echo "Installing build dependencies..."
    if command -v mk-build-deps >/dev/null 2>&1; then
        # Use mk-build-deps (creates and installs a metapackage)
        sudo mk-build-deps --install --remove --tool='apt-get -o Debug::pkgProblemResolver=yes --no-install-recommends --yes' debian/control
    elif command -v apt-get >/dev/null 2>&1; then
        # Fallback: use apt-get build-dep
        echo "mk-build-deps not found, trying apt-get build-dep..."
        sudo apt-get update
        sudo apt-get build-dep -y .
    else
        echo "Warning: Cannot install build dependencies automatically"
        echo "Install manually: sudo apt-get build-dep ."
    fi
fi

echo "Building package..."
debuild -us -uc

echo ""
echo "Build complete! Check ../ for .deb files"
