#!/bin/bash
# Auto-detect latest wsclean upstream version and fetch with submodules
# Usage: ./debian/fetch-latest.sh [--distro noble] [--filter v3]
# Example: ./debian/fetch-latest.sh                    # Gets latest
# Example: ./debian/fetch-latest.sh --filter v3        # Gets latest v3.x
# Example: ./debian/fetch-latest.sh --distro focal     # Latest for focal

set -e

UPSTREAM_URL="https://gitlab.com/aroffringa/wsclean.git"
DISTRO="${DISTRO:-noble}"
VERSION_FILTER="${1:-}"
WORK_DIR=""

echo "WSClean Auto-Fetch Latest Version"
echo ""

# Function to parse version from tag
get_latest_version() {
    local filter="$1"

    echo "Querying GitLab for latest tags..." >&2
    
    # Fetch all tags from upstream, sort by version
    # Filter format: "v3" matches v3.0, v3.1, etc
    git ls-remote --tags "$UPSTREAM_URL" | \
        grep "refs/tags/v" | \
        awk -F'/' '{print $NF}' | \
        sed 's/\^{}//' | \
        sort -V | \
        tail -10 > /tmp/wsclean_tags.txt
    
    if [ -z "$filter" ]; then
        # Get absolute latest
        LATEST=$(tail -1 /tmp/wsclean_tags.txt)
    else
        # Get latest matching filter (e.g., "v3" returns latest v3.x)
        LATEST=$(grep "^${filter}" /tmp/wsclean_tags.txt | tail -1)
    fi
    
    if [ -z "$LATEST" ]; then
        echo "No matching tags found" >&2
        echo "Available tags:" >&2
        cat /tmp/wsclean_tags.txt >&2
        exit 1
    fi
    
    # Remove 'v' prefix for VERSION variable
    echo "$LATEST"
}

# Get latest version
LATEST_TAG=$(get_latest_version "$VERSION_FILTER")
VERSION="${LATEST_TAG#v}"  # Strip 'v' prefix

echo "Latest version: $VERSION (tag: $LATEST_TAG)"
echo ""

PKG_NAME="wsclean"
WORK_DIR="${PKG_NAME}-${VERSION}"
TARBALL="${PKG_NAME}_${VERSION}.orig.tar.gz"

# Clean previous builds
if [ -d "${WORK_DIR}" ]; then
    echo "Removing previous build directory..."
    rm -rf "${WORK_DIR}"
fi

if [ -f "${TARBALL}" ]; then
    echo "Removing previous tarball..."
    rm -f "${TARBALL}"
fi

# Clone with submodules
echo "Cloning upstream ${LATEST_TAG} with recursive submodules..."
git clone \
    --recursive \
    --depth 1 \
    --branch "${LATEST_TAG}" \
    "$UPSTREAM_URL" "${WORK_DIR}"

cd "${WORK_DIR}"

echo "Verifying submodule initialization..."
git submodule update --init --recursive

# Remove Git metadata
echo "Stripping Git directories..."
find . -name ".git*" -type d -print0 | xargs -0 rm -rf
find . -name ".gitignore" -type f -delete
find . -name ".gitattributes" -type f -delete
find . -name ".gitlab-ci.yml" -type f -delete

# Create tarball
cd ..
echo "Creating orig.tar.gz..."
tar \
    --exclude='*.swp' \
    --exclude='*.swo' \
    --exclude='.DS_Store' \
    -czf "${TARBALL}" "${WORK_DIR}"

# Verify
TARBALL_SIZE=$(du -h "${TARBALL}" | cut -f1)
TARBALL_FILES=$(tar -tzf "${TARBALL}" | wc -l)

echo ""
echo "SUCCESS! Auto-fetched and prepared:"
echo "   Version:  ${VERSION}"
echo "   Tag:      ${LATEST_TAG}"
echo "   File:     ${TARBALL}"
echo "   Size:     ${TARBALL_SIZE}"
echo "   Files:    ${TARBALL_FILES}"
echo ""
echo "Next step:"
echo "   $ debian/build-and-sign.sh ${VERSION} ${DISTRO}"
echo ""
