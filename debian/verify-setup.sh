#!/bin/bash
# Verification script to test the full auto-build workflow
# This script tests everything WITHOUT actually building (dry-run)

set -e

echo ""
echo "================================================================"
echo "        WSClean Auto-Build Verification (Dry-Run)"
echo "   This tests the workflow without consuming bandwidth"
echo "================================================================"
echo ""

echo "Verification Checklist:"
echo ""

# Check 1: Scripts exist and are executable
echo "[*] Checking scripts..."
for script in fetch-latest.sh build-and-sign.sh full-build.sh; do
    if [ -x "debian/$script" ]; then
        echo "  [OK] debian/$script (executable)"
    else
        echo "  [FAIL] debian/$script (missing or not executable!)"
        exit 1
    fi
done
echo ""

# Check 2: Can connect to upstream
echo "[*] Testing upstream access..."
if git ls-remote --tags https://gitlab.com/aroffringa/wsclean.git >/dev/null 2>&1; then
    echo "  [OK] Can reach GitLab upstream"
else
    echo "  [FAIL] Cannot reach GitLab!"
    exit 1
fi
echo ""

# Check 3: Get latest version
echo "[*] Auto-detecting latest version..."
LATEST_TAG=$(git ls-remote --tags https://gitlab.com/aroffringa/wsclean.git | \
    grep "refs/tags/v" | \
    awk -F'/' '{print $NF}' | \
    sed 's/\^{}//' | \
    sort -V | \
    tail -1)

VERSION="${LATEST_TAG#v}"
echo "  [OK] Latest version: ${VERSION} (tag: ${LATEST_TAG})"
echo ""

# Check 4: Required tools
echo "[*] Checking required tools..."
for tool in git cmake tar gzip; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "  [OK] $tool"
    else
        echo "  [FAIL] $tool not found!"
        exit 1
    fi
done
echo ""

# Check 5: Debian tools
echo "[*] Checking Debian packaging tools..."
missing=0
for tool in debuild dch debsign; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "  [OK] $tool"
    else
        echo "  [FAIL] $tool (missing - install: sudo apt install devscripts)"
        missing=$((missing+1))
    fi
done

if [ "$missing" -gt 0 ]; then
    echo ""
    echo "Missing Debian tools! Install with:"
    echo "   sudo apt install build-essential debhelper devscripts dh-make"
    exit 1
fi
echo ""

# Check 6: GPG setup
echo "[*] Checking GPG signing setup..."
if command -v gpg >/dev/null 2>&1; then
    KEY_COUNT=$(gpg --list-secret-keys 2>/dev/null | grep "sec   " | wc -l)
    if [ "$KEY_COUNT" -gt 0 ] 2>/dev/null; then
        echo "  [OK] GPG ($KEY_COUNT signing key(s) available)"
        gpg --list-secret-keys --keyid-format short 2>/dev/null | grep uid | head -1 | sed 's/^/    /'
    else
        echo "  [INFO] No GPG keys found (signing will require manual input)"
    fi
else
    echo "  [FAIL] GPG not installed"
    exit 1
fi
echo ""

# Check 7: Disk space
echo "[*] Checking free space..."
free_space=$(df /workspaces 2>/dev/null | awk '{if(NR==2) print int($4)}')
if [ -n "$free_space" ] && [ "$free_space" -gt 1000000 ]; then  # > 1GB
    free_gb=$((free_space / 1048576))
    echo "  [OK] ${free_gb}GB free (sufficient)"
else
    echo "  [WARN] Low disk space or cannot determine!"
fi
echo ""

# Check 8: Network (can clone from upstream)
echo "[*] Testing clone capability..."
if timeout 5 git clone --depth 1 --branch v3.4 --no-checkout \
    https://gitlab.com/aroffringa/wsclean.git /tmp/wsclean_test 2>/dev/null; then
    rm -rf /tmp/wsclean_test
    echo "  [OK] Can clone from upstream"
else
    echo "  [FAIL] Clone test failed (network issue?)"
    exit 1
fi
echo ""

# Final summary
echo "================================================================"
echo "              ALL CHECKS PASSED!"
echo "================================================================"
echo ""
echo "Ready to build!"
echo ""
echo "Next step: Run the full workflow"
echo ""
echo "  Option 1: Auto-fetch latest (recommended)"
echo "  $ debian/full-build.sh --distro noble"
echo ""
echo "  Option 2: Fetch only (requires manual build)"
echo "  $ debian/fetch-latest.sh"
echo ""
echo "  Option 3: Build specific version"
echo "  $ debian/fetch-upstream.sh 3.6"
echo "  $ debian/build-and-sign.sh 3.6 noble"
echo ""
echo "Current upstream status:"
echo "  Latest:  $VERSION (auto-detected)"
echo "  Release: $(echo $LATEST_TAG | tr -d 'v')"
echo ""
