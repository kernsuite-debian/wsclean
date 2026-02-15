#!/bin/bash
# Build, sign, and prepare local upload for wsclean debian package
# Usage: ./debian/build-and-sign.sh <VERSION> [DISTRO] [GPG_KEY]
# Example: ./debian/build-and-sign.sh 3.6 noble
# Example: ./debian/build-and-sign.sh 3.6 noble ABC123DEF

set -e

VERSION="${1}"
DISTRO="${2:-noble}"
GPG_KEY="${3}"

if [ -z "$VERSION" ]; then
    echo "Usage: $0 <VERSION> [DISTRO] [GPG_KEY_ID]"
    echo ""
    echo "Example:"
    echo "  $0 3.6 noble                    # Build, sign with default key"
    echo "  $0 3.6 noble 12345678          # Build, sign with specific GPG key"
    exit 1
fi

PKG_NAME="wsclean"
TARBALL="${PKG_NAME}_${VERSION}.orig.tar.gz"
BUILD_DIR="../build-area"
UPLOAD_DIR="../wsclean-debs-${VERSION}"

echo ""
echo "================================================================"
echo "     WSClean Build, Sign & Local Upload Workflow"
echo "================================================================"
echo ""

echo "Configuration:"
echo "   Version:    ${VERSION}"
echo "   Distro:     ${DISTRO}"
echo "   Tarball:    ${TARBALL}"
echo "   Build dir:  ${BUILD_DIR}"
echo "   Upload dir: ${UPLOAD_DIR}"
if [ -n "$GPG_KEY" ]; then
    echo "   GPG Key:    ${GPG_KEY}"
fi
echo ""

# Verify tarball exists
if [ ! -f "../${TARBALL}" ]; then
    echo "ERROR: Tarball not found: ../${TARBALL}"
    echo ""
    echo "First, run: debian/fetch-latest.sh"
    exit 1
fi

# Step 1: Update changelog
echo "----------------------------------------------------------------"
echo "Step 1: Updating debian/changelog..."
echo "----------------------------------------------------------------"

if command -v dch &>/dev/null; then
    # Check if version already exists in changelog
    CURRENT_TOP=$(head -1 debian/changelog | grep -oP '\(\K[^)]+')
    NEW_VERSION=""
    
    if echo "$CURRENT_TOP" | grep -q "^${VERSION}-"; then
        # Version exists, increment kern number
        if echo "$CURRENT_TOP" | grep -qP "${VERSION}-\d+kern\d+"; then
            # Extract and increment kern number
            DEBIAN_REV=$(echo "$CURRENT_TOP" | grep -oP "${VERSION}-\K\d+")
            KERN_NUM=$(echo "$CURRENT_TOP" | grep -oP "kern\K\d+")
            NEW_KERN=$((KERN_NUM + 1))
            NEW_VERSION="${VERSION}-${DEBIAN_REV}kern${NEW_KERN}"
            echo "[INFO] Version ${VERSION} exists, incrementing to ${NEW_VERSION}"
        else
            NEW_VERSION="${VERSION}-1kern1"
        fi
    else
        # New upstream version
        NEW_VERSION="${VERSION}-1kern1"
        echo "[INFO] New upstream version ${VERSION}"
    fi
    
    DEBFULLNAME="KERN packaging" DEBEMAIL="packaging@kernsuite.info" \
    dch --newversion "${NEW_VERSION}" \
        --distribution "${DISTRO}" \
        "New upstream release v${VERSION}"
    echo "[OK] Changelog updated to ${NEW_VERSION}"
else
    echo "[WARN] dch not available, skipping automatic changelog update"
    echo "   Manual edit debian/changelog with:"
    echo "   ${PKG_NAME} (${VERSION}-1kern1) ${DISTRO}; urgency=medium"
fi
echo ""

# Step 2: Import into git-buildpackage (optional but recommended)
echo "----------------------------------------------------------------"
echo "Step 2: Building package..."
echo "----------------------------------------------------------------"

# Create build directory
mkdir -p "${BUILD_DIR}"

# Build (unsigned first for testing)
echo "Building (unsigned test build)..."
debuild --no-tgz-check -us -uc 2>&1 | tail -30

# Check if build succeeded
if [ ! -f "../${PKG_NAME}_${VERSION}-1kern1_"*.deb ]; then
    echo ""
    echo "[FAIL] Build failed! Check output above for errors."
    exit 1
fi

echo ""
echo "[OK] Build succeeded!"
echo ""

# Step 3: List generated files
echo "----------------------------------------------------------------"
echo "Step 3: Generated files (ready for signing)..."
echo "----------------------------------------------------------------"

ls -lh ../${PKG_NAME}_${VERSION}-1kern1* | grep -v ".build"
echo ""

# Step 4: Sign files
echo "----------------------------------------------------------------"
echo "Step 4: Signing files..."
echo "----------------------------------------------------------------"

# Sign the .changes file (this signs all referenced files)
CHANGES_FILE="../${PKG_NAME}_${VERSION}-1kern1_"*".changes"

if [ -f "$CHANGES_FILE" ]; then
    # Check if already signed
    if grep -q "^BEGIN PGP SIGNATURE" "$CHANGES_FILE"; then
        echo "[WARN] Changes file already signed"
    else
        echo "Signing ${CHANGES_FILE##*/}..."
        if [ -n "$GPG_KEY" ]; then
            debsign -k "$GPG_KEY" "$CHANGES_FILE"
        else
            debsign "$CHANGES_FILE"
        fi
        echo "[OK] Files signed"
    fi
else
    echo "[FAIL] Changes file not found!"
    exit 1
fi
echo ""

# Step 5: Create upload directory
echo "----------------------------------------------------------------"
echo "Step 5: Preparing local upload directory..."
echo "----------------------------------------------------------------"

mkdir -p "${UPLOAD_DIR}"

# Copy all relevant files
cp ../${PKG_NAME}_${VERSION}-1kern1_*.deb "${UPLOAD_DIR}/" 2>/dev/null || true
cp ../${PKG_NAME}_${VERSION}-1kern1_*.changes "${UPLOAD_DIR}/"
cp ../${PKG_NAME}_${VERSION}-1kern1_*.dsc "${UPLOAD_DIR}/"
cp ../${TARBALL} "${UPLOAD_DIR}/"

echo "Files prepared in: ${UPLOAD_DIR}"
ls -lh "${UPLOAD_DIR}/"
echo ""

# Step 6: Run lintian checks
echo "----------------------------------------------------------------"
echo "Step 6: Running quality checks (lintian)..."
echo "----------------------------------------------------------------"

if command -v lintian &>/dev/null; then
    echo "Checking .deb packages..."
    lintian "${UPLOAD_DIR}"/${PKG_NAME}_${VERSION}-1kern1_*.deb 2>&1 | head -20 || true
    
    echo ""
    echo "Checking .dsc source..."
    lintian "${UPLOAD_DIR}"/${PKG_NAME}_${VERSION}-1kern1_*.dsc 2>&1 | head -20 || true
else
    echo "[WARN] lintian not available, skipping package checks"
fi
echo ""

# Step 7: Verification
echo "----------------------------------------------------------------"
echo "Step 7: Verification Summary"
echo "----------------------------------------------------------------"

echo "[OK] Package built: $(ls -1 ${UPLOAD_DIR}/${PKG_NAME}_${VERSION}-1kern1_*.deb | wc -l) .deb file(s)"
echo "[OK] Signed: $(grep -l 'BEGIN PGP SIGNATURE' ${UPLOAD_DIR}/*.changes)"
echo "[OK] Ready for upload: ${UPLOAD_DIR}"
echo ""

# Step 8: Show upload options
echo "----------------------------------------------------------------"
echo "Step 8: Next Steps - Upload Options"
echo "----------------------------------------------------------------"

CHANGES_FILE_SHORT="${PKG_NAME}_${VERSION}-1kern1_"*.changes

echo ""
echo "Option A: Upload to local apt repository"
echo "  $ cp ${UPLOAD_DIR}/* /path/to/local/apt/repo/pool/main/${PKG_NAME}/"
echo "  $ cd /path/to/local/apt/repo && reprepro -Vb . includedeb ${DISTRO} ../${PKG_NAME}*.deb"
echo ""

echo "Option B: Upload to remote (requires SSH access)"
echo "  $ dput -u ${UPLOAD_DIR}/${CHANGES_FILE_SHORT}"
echo ""

echo "Option C: Test locally before upload"
echo "  $ sudo dpkg -i ${UPLOAD_DIR}/${PKG_NAME}_${VERSION}-1kern1_*.deb"
echo "  $ wsclean --version"
echo ""

echo "Option D: Create an archive for distribution"
echo "  $ tar -czf wsclean-${VERSION}-debs.tar.gz ${UPLOAD_DIR}/"
echo ""

echo "Ready to upload!"
echo ""
