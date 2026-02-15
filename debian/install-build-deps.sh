#!/bin/bash
# Install build dependencies for wsclean
# Run this before building locally

set -e

echo "Installing build dependencies for wsclean..."

# Check if mk-build-deps is available
if ! command -v mk-build-deps >/dev/null 2>&1; then
    echo "Installing devscripts (provides mk-build-deps)..."
    sudo apt-get update
    sudo apt-get install -y devscripts equivs
fi

# Use mk-build-deps to install all dependencies
echo "Using mk-build-deps to install dependencies..."
cd "$(dirname "$0")/.."
sudo mk-build-deps --install --remove \
    --tool='apt-get -o Debug::pkgProblemResolver=yes --no-install-recommends --yes' \
    debian/control

echo ""
echo "Build dependencies installed successfully!"
echo ""
echo "You can now build with:"
echo "  debuild -us -uc"
echo ""
