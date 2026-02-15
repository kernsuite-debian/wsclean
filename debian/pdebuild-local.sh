#!/bin/bash
# Prepare and run pdebuild - for use after fresh clone
# Usage: ./debian/pdebuild-local.sh [VERSION]

set -e

# Auto-detect version from changelog if not provided
if [ -z "$1" ]; then
    VERSION=$(dpkg-parsechangelog -S Version | cut -d'-' -f1)
    echo "Auto-detected version: $VERSION"
else
    VERSION="$1"
fi

PKG_NAME="wsclean"
TARBALL="../${PKG_NAME}_${VERSION}.orig.tar.gz"

echo "================================================================"
echo "Prepare and Run pdebuild for ${PKG_NAME} ${VERSION}"
echo "================================================================"
echo ""

# Check if tarball exists
if [ ! -f "$TARBALL" ]; then
    echo "Upstream tarball not found: $TARBALL"
    echo ""
    echo "Fetching upstream version ${VERSION}..."
    debian/fetch-upstream.sh "$VERSION"
    echo ""
fi

# Verify tarball exists now
if [ ! -f "$TARBALL" ]; then
    echo "Error: Failed to create tarball!"
    exit 1
fi

echo "Tarball found: $TARBALL"
echo "Running pdebuild..."
echo ""

# Run pdebuild
pdebuild "$@"

echo ""
echo "================================================================"
echo "Build complete!"
echo "================================================================"
echo ""
echo "Results in: /var/cache/pbuilder/result/"
echo ""
