#!/bin/bash
# Test installation of built .deb packages
# This removes build dependencies first to verify runtime deps are correct

set -e

VERSION="${1}"
if [ -z "$VERSION" ]; then
    echo "Usage: $0 <VERSION>"
    echo "Example: $0 3.6"
    exit 1
fi

echo "================================================================"
echo "Test Installation: Verify Runtime Dependencies"
echo "================================================================"
echo ""

# Remove build dependencies (mark as auto-installed, then autoremove)
echo "Step 1: Cleaning up build dependencies..."
if dpkg -l | grep -q "wsclean-build-deps"; then
    sudo apt-get remove -y wsclean-build-deps
fi
sudo apt-get autoremove -y

echo ""
echo "Step 2: Installing built packages..."
DEB_FILES=$(ls ../wsclean*_${VERSION}*.deb 2>/dev/null | grep -v build-deps || true)

if [ -z "$DEB_FILES" ]; then
    echo "Error: No .deb files found for version ${VERSION}"
    echo "Build the package first with: debian/build-from-tarball.sh ${VERSION}"
    exit 1
fi

echo "Found packages:"
echo "$DEB_FILES" | sed 's/^/  /'
echo ""

# Try to install
sudo dpkg -i $DEB_FILES || {
    echo ""
    echo "Missing runtime dependencies detected!"
    echo "Attempting to fix with apt-get..."
    sudo apt-get install -f -y
}

echo ""
echo "Step 3: Testing installed package..."
if command -v wsclean >/dev/null 2>&1; then
    echo "[OK] wsclean is installed"
    wsclean --version
    echo ""
    echo "[SUCCESS] Installation test passed!"
else
    echo "[FAIL] wsclean not found in PATH"
    exit 1
fi

echo ""
echo "To uninstall:"
echo "  sudo apt-get remove wsclean wsclean-dev libwsclean2"
echo ""
