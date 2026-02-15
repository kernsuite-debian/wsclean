#!/bin/bash
# Complete local workflow: Clone → Test → Sign → Upload
# Run this on your LOCAL machine after pushing changes to GitHub
#
# Usage: 
#   ./debian/local-build-and-upload.sh [PPA1] [PPA2...]
# Example:
#   ./debian/local-build-and-upload.sh ppa:kernsuite/kern-dev ppa:kernsuite/kern-10

set -e

# Default to both kern-dev and kern-10 if no PPA specified
if [ $# -eq 0 ]; then
    PPAS=("ppa:kernsuite/kern-dev" "ppa:kernsuite/kern-10")
else
    PPAS=("$@")
fi

VERSION=$(dpkg-parsechangelog -S Version | cut -d'-' -f1)
FULL_VERSION=$(dpkg-parsechangelog -S Version)
PKG_NAME="wsclean"

echo "================================================================"
echo "WSClean Local Build → Sign → PPA Upload Workflow"
echo "================================================================"
echo ""
echo "Package:  ${PKG_NAME}"
echo "Version:  ${FULL_VERSION}"
echo "PPAs:     ${PPAS[@]}"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
fi

# Step 1: Check for upstream tarball
echo ""
echo "----------------------------------------------------------------"
echo "Step 1: Verify/Create Upstream Tarball"
echo "----------------------------------------------------------------"
TARBALL="../${PKG_NAME}_${VERSION}.orig.tar.gz"
if [ ! -f "$TARBALL" ]; then
    echo "Tarball not found, fetching upstream..."
    debian/fetch-upstream.sh "$VERSION"
else
    echo "[OK] Tarball exists: $TARBALL"
fi

# Step 2: Clean room test with pdebuild (optional but recommended)
echo ""
echo "----------------------------------------------------------------"
echo "Step 2: Clean Room Build Test (pdebuild)"
echo "----------------------------------------------------------------"
read -p "Run pdebuild test? (recommended) (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Running pdebuild..."
    if pdebuild; then
        echo "[OK] pdebuild succeeded"
        echo "Results in: /var/cache/pbuilder/result/"
    else
        echo "[FAIL] pdebuild failed!"
        read -p "Continue anyway? (y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi
else
    echo "Skipping pdebuild..."
fi

# Step 3: Build source package
echo ""
echo "----------------------------------------------------------------"
echo "Step 3: Build Source Package for PPA"
echo "----------------------------------------------------------------"
echo "Building source package..."
debuild -S -sa

# Check if source package was created
CHANGES_FILE="../${PKG_NAME}_${FULL_VERSION}_source.changes"
if [ ! -f "$CHANGES_FILE" ]; then
    echo "[FAIL] Source package not created!"
    exit 1
fi

echo "[OK] Source package created: ${CHANGES_FILE}"

# Step 4: Verify signing
echo ""
echo "----------------------------------------------------------------"
echo "Step 4: Verify Package is Signed"
echo "----------------------------------------------------------------"
if grep -q "BEGIN PGP SIGNATURE" "$CHANGES_FILE"; then
    echo "[OK] Package is signed"
    echo ""
    echo "Signed by:"
    gpg --verify "$CHANGES_FILE" 2>&1 | grep "Good signature" || echo "  (signature details above)"
else
    echo "[FAIL] Package is NOT signed!"
    echo "You need a GPG key configured."
    echo ""
    echo "To sign manually:"
    echo "  debsign ${CHANGES_FILE}"
    exit 1
fi

# Step 5: Upload to PPA
echo ""
echo "----------------------------------------------------------------"
echo "Step 5: Upload to PPA(s)"
echo "----------------------------------------------------------------"
echo "About to upload to:"
for ppa in "${PPAS[@]}"; do
    echo "  - $ppa"
done
echo ""
echo "Files:"
ls -lh ../${PKG_NAME}_${FULL_VERSION}* | grep -v ".build"
echo ""
read -p "Upload now? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    for ppa in "${PPAS[@]}"; do
        echo ""
        echo "Uploading to $ppa..."
        dput "$ppa" "$CHANGES_FILE"
        echo "[OK] Uploaded to $ppa"
    done
    echo ""
    echo "[SUCCESS] Uploaded to all PPAs"
    echo ""
    echo "Track build status at:"
    echo "  https://launchpad.net/~kernsuite/+archive/ubuntu/kern-dev/+packages"
    echo "  https://launchpad.net/~kernsuite/+archive/ubuntu/kern-10/+packages"
else
    echo "Upload skipped."
    echo ""
    echo "To upload manually:"
    for ppa in "${PPAS[@]}"; do
        echo "  dput $ppa ${CHANGES_FILE}"
    done
fi

echo ""
echo "================================================================"
echo "                      Workflow Complete!"
echo "================================================================"
echo ""
