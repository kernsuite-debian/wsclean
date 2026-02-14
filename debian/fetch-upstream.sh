#!/bin/bash
# AI-Ready script to handle wsclean submodules for Debian packaging
# Usage: ./fetch-upstream.sh 3.5
# This script automates the clone-submodule-tarball workflow

set -e

VERSION="${1}"
if [ -z "$VERSION" ]; then
    echo "Usage: $0 <VERSION>"
    echo "Example: $0 3.5"
    exit 1
fi

UPSTREAM_URL="https://gitlab.com/aroffringa/wsclean.git"
PKG_NAME="wsclean"
WORK_DIR="${PKG_NAME}-${VERSION}"
TARBALL="${PKG_NAME}_${VERSION}.orig.tar.gz"

echo " Starting wsclean Debian upstream fetch workflow..."
echo "   Version: ${VERSION}"
echo "   Target: ${TARBALL}"

# Step 1: Clean previous builds
if [ -d "${WORK_DIR}" ]; then
    echo "  Removing previous build directory..."
    rm -rf "${WORK_DIR}"
fi

if [ -f "${TARBALL}" ]; then
    echo "  Removing previous tarball..."
    rm -f "${TARBALL}"
fi

# Step 2: Clone with submodules
echo "Cloning upstream repository with submodules..."
git clone \
    --recursive \
    --depth 1 \
    --branch "v${VERSION}" \
    "$UPSTREAM_URL" "${WORK_DIR}"

cd "${WORK_DIR}"

# Step 3: Ensure all submodules are initialized (safety check)
echo "Verifying submodule initialization..."
git submodule update --init --recursive

# Step 4: Remove Git metadata
echo "Stripping Git directories and metadata..."
find . -name ".git*" -type d -print0 | xargs -0 rm -rf
find . -name ".gitignore" -type f -delete
find . -name ".gitattributes" -type f -delete
find . -name ".gitlab-ci.yml" -type f -delete

# Step 5: Create tarball
cd ..
echo "Creating orig tarball..."
tar \
    --exclude='*.swp' \
    --exclude='*.swo' \
    --exclude='.DS_Store' \
    -czf "${TARBALL}" "${WORK_DIR}"

# Step 6: Verify integrity
TARBALL_SIZE=$(du -h "${TARBALL}" | cut -f1)
TARBALL_FILES=$(tar -tzf "${TARBALL}" | wc -l)

echo ""
echo "✅ SUCCESS! Tarball created:"
echo "   File: ${TARBALL}"
echo "   Size: ${TARBALL_SIZE}"
echo "   Files: ${TARBALL_FILES}"
echo ""
echo "Next steps:"
echo "   1. Update debian/changelog with new version"
echo "   2. Run: debuild -us -uc"
echo "   3. Check with: lintian ../wsclean_${VERSION}-*.deb"
echo ""
