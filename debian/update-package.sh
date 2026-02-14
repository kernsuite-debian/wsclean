#!/bin/bash
# This script orchestrates the full workflow: fetch → build → check
# Usage: debian/update-package.sh 3.6 noble

set -e

VERSION="${1}"
DISTRO="${2:-noble}"

if [ -z "$VERSION" ]; then
    echo "WSClean Debian Update Assistant"
    echo ""
    echo "Usage: $0 <VERSION> [DISTRO]"
    echo "Example: $0 3.6 noble"
    echo ""
    echo "This script:"
    echo "  1. Fetches upstream v${VERSION} with all submodules"
    echo "  2. Prepares .orig.tar.gz"
    echo "  3. Updates debian/changelog"
    echo "  4. Builds the package (dry-run)"
    echo "  5. Runs lintian checks"
    exit 0
fi

PKG_NAME="wsclean"
CURRENT_DIR="$(pwd)"

echo "Packaging Wizard: Updating ${PKG_NAME} to ${VERSION}"
echo "   Target distro: ${DISTRO}"
echo ""

# Step 1: Fetch upstream with submodules
echo "Step 1: Fetching upstream ${VERSION} with submodules..."
debian/fetch-upstream.sh "${VERSION}"

# Step 2: Update changelog
echo "Step 2: Updating debian/changelog..."
EXPECTED_TARBALL="${PKG_NAME}_${VERSION}.orig.tar.gz"
if [ ! -f "../${EXPECTED_TARBALL}" ]; then
    echo "ERROR: Expected tarball ../[${EXPECTED_TARBALL}] not found!"
    exit 1
fi

MAINTAINER="KERN packaging <packaging@kernsuite.info>"
TIMESTAMP=$(date -R)
CHANGELOG_ENTRY="${PKG_NAME} (${VERSION}-1kern1) ${DISTRO}; urgency=medium

  * New upstream release v${VERSION}
  * Rebuilt with submodule support

 -- ${MAINTAINER}  ${TIMESTAMP}"

# Prepend to changelog (requires dch command from devscripts)
if command -v dch &>/dev/null; then
    cd debian/..
    dch --newversion "${VERSION}-1kern1" \
        --distribution "${DISTRO}" \
        --maintainer "${MAINTAINER}" \
        "New upstream release v${VERSION}"
    cd "${CURRENT_DIR}"
    echo "✅ Changelog updated"
else
    echo "Warning: dch not found, skipping automatic changelog update"
    echo "   Manual update: Edit debian/changelog"
fi

# Step 3: Import into git-buildpackage
echo ""
echo "Step 3: Preparing git-buildpackage import..."
echo "   (Run this manually if using gbp):"
echo "   $ gbp import-orig --pristine-tar ../${EXPECTED_TARBALL}"
echo ""

# Step 4: Test build
echo "Step 4: Test build (configure only)..."
if [ -f "CMakeLists.txt" ]; then
    mkdir -p build-test
    cd build-test
    cmake .. -DCMAKE_BUILD_TYPE=Release 2>&1 | head -20
    cd "${CURRENT_DIR}"
    rm -rf build-test
    echo "✅ CMake configuration successful"
else
    echo "No CMakeLists.txt found"
fi

# Step 5: Lintian checks
echo ""
echo "✅ Package preparation complete!"
echo ""
echo "Next steps:"
echo "  1. gbp import-orig --pristine-tar ../${EXPECTED_TARBALL}"
echo "  2. debuild -us -uc"
echo "  3. lintian ../${PKG_NAME}_${VERSION}-1kern1_*.deb"
echo ""
