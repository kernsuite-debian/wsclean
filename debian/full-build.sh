#!/bin/bash
# Full workflow: Auto-fetch latest + Build + Sign + Prepare for upload
# Usage: ./debian/full-build.sh [--distro noble] [--sign GPG_KEY]
#
# This is the all-in-one command that requires NO version number

set -e

DISTRO="noble"
GPG_KEY=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --distro) DISTRO="$2"; shift 2 ;;
        --sign) GPG_KEY="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  WSClean Full Workflow: Fetch Latest → Build → Sign → Upload║"
echo "║                                                            ║"
echo "║  This script auto-detects the latest version and builds!  ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Step 1: Fetch latest
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1: Auto-Detect and Fetch Latest Version"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

debian/fetch-latest.sh

# Extract version from generated tarball
TARBALL=$(ls -1 wsclean_*.orig.tar.gz | tail -1)
VERSION="${TARBALL#wsclean_}"
VERSION="${VERSION%.orig.tar.gz}"

echo "🎯 Detected version: ${VERSION}"
echo ""

# Step 2: Build and sign
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 2: Build, Sign, and Prepare Upload"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -n "$GPG_KEY" ]; then
    debian/build-and-sign.sh "$VERSION" "$DISTRO" "$GPG_KEY"
else
    debian/build-and-sign.sh "$VERSION" "$DISTRO"
fi

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║                    🎉 ALL DONE! 🎉                         ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "✨ Summary:"
echo "  Version:      $VERSION"
echo "  Distribution: $DISTRO"
echo "  Upload dir:   ../wsclean-debs-${VERSION}"
echo ""
echo "Next: Check ../wsclean-debs-${VERSION} for .deb files"
echo ""
