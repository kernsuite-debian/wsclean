#!/bin/bash
# Debian repack script: Handles submodules when uscan downloads upstream
# This script is called by uscan or gbp when a new upstream tarball is found
# It extracts, adds submodules, and repackages as .orig.tar.gz

set -e

VERSION="${1}"
TARBALL="${2:-wsclean_${VERSION}.tar.gz}"

if [ -z "$VERSION" ]; then
    echo "Usage: $0 <VERSION> [TARBALL]"
    echo ""
    echo "This script repackages a uscan-downloaded tarball to include git submodules."
    echo "Call from gbp with: gbp import-orig --uscan --pristine-tar"
    exit 1
fi

PKG_NAME="wsclean"
UPSTREAM_URL="https://gitlab.com/aroffringa/wsclean.git"
WORK_DIR="${PKG_NAME}-${VERSION}"

echo "Repack: Enhancing uscan tarball with git submodules..."
echo "   Version: ${VERSION}"
echo "   Input:   ${TARBALL}"

# Step 1: Extract the uscan-provided tarball
if [ ! -f "${TARBALL}" ]; then
    echo "ERROR: Tarball not found: ${TARBALL}"
    exit 1
fi

echo "Extracting base tarball..."
mkdir -p "${WORK_DIR}"
tar -xzf "${TARBALL}" -C "${WORK_DIR}" --strip-components=1

# Step 2: Clone fresh repo and copy submodules
echo "Fetching git submodules..."
TEMP_REPO=$(mktemp -d)
git clone --recursive --branch "v${VERSION}" "$UPSTREAM_URL" "$TEMP_REPO"

# Copy .gitmodules if present (to preserve submodule information)
if [ -f "$TEMP_REPO/.gitmodules" ]; then
    cp "$TEMP_REPO/.gitmodules" "${WORK_DIR}/"
fi

# Copy all submodule directories
for submodule_dir in "$TEMP_REPO"/external/*; do
    if [ -d "$submodule_dir" ]; then
        submodule_name=$(basename "$submodule_dir")
        echo "  • Copying submodule: $submodule_name"
        if [ -d "${WORK_DIR}/external/$submodule_name" ]; then
            rm -rf "${WORK_DIR}/external/$submodule_name"
        fi
        cp -r "$submodule_dir" "${WORK_DIR}/external/"
    fi
done

rm -rf "$TEMP_REPO"

# Step 3: Clean git metadata
echo "Removing git metadata..."
find "${WORK_DIR}" -name ".git*" -type d -print0 | xargs -0 rm -rf
find "${WORK_DIR}" -name ".gitignore" -delete
find "${WORK_DIR}" -name ".gitlab-ci.yml" -delete

# Step 4: Repackage as .orig.tar.gz
ORIG_TARBALL="${PKG_NAME}_${VERSION}.orig.tar.gz"
echo "📦 Creating orig.tar.gz..."
tar -czf "${ORIG_TARBALL}" "${WORK_DIR}"

# Clean build dir
rm -rf "${WORK_DIR}"

echo "✅ SUCCESS!"
echo "   Output: ${ORIG_TARBALL}"
echo ""
echo "Next steps:"
echo "   1. gbp import-orig --pristine-tar ${ORIG_TARBALL}"
echo "   2. Update debian/changelog"
echo "   3. debuild -us -uc"
