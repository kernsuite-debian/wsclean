# WSClean Debian Packaging Guide

## Overview

This guide documents the packaging workflow for wsclean. Instead of manual clone-rename-delete cycles, we use automated scripts that handle submodules correctly.

**Key components:**
- `debian/fetch-upstream.sh` — Fetches upstream with submodules, creates .orig.tar.gz
- `debian/repack.sh` — Enhances uscan tarballs with submodules (for traditional watch file)
- `debian/update-package.sh` — Full orchestration: fetch + build + check
- `.devcontainer/` — Pre-configured Codespace environment

---

## Quick Start (First Time)

### 1. Set Up Codespace Environment

```bash
# In your Codespace, the environment is auto-configured when you open the repo
# Verify tools are installed:
which debuild dch lintian git-buildpackage

# If running locally, install dependencies:
sudo apt install build-essential debhelper devscripts dh-make \
    git-buildpackage pristine-tar pkg-config cmake doxygen
```

### 2. Prepare Upstream with Submodules

```bash
cd /path/to/wsclean-repo

# Fetch version 3.6 with all submodules
debian/fetch-upstream.sh 3.6

# This creates: wsclean_3.6.orig.tar.gz in parent directory
```

### 3. Update Package Metadata

```bash
# One-command full workflow (fetch + changelog + build check)
debian/update-package.sh 3.6 noble

# Manual steps if you prefer:
dch --newversion 3.6-1kern1 \
    --distribution noble \
    "New upstream release v3.6"
```

### 4. Build & Check

```bash
# Build unsigned (for testing)
debuild -us -uc

# Run lintian checks
lintian ../wsclean_3.6-1kern1_*.deb
```

---

## For AI Agents: Mission Commands

### Command 1: "Upgrade Package"

**Prompt for AI (Copilot/Agent):**
```
I need to update wsclean to version 3.6 in Debian format.

1. Run: debian/update-package.sh 3.6 noble
2. Check the output for any build warnings
3. If CMake succeeds, tell me the next commands to build the .deb
4. If there are errors, suggest fixes based on the CMakeLists.txt
```

**What the AI will do:**
- Execute fetch-upstream.sh (submodules included)
- Update debian/changelog via dch
- Do a test CMake configure
- Provide lintian instructions

---

### Command 2: "Check Build Dependencies"

**Prompt for AI:**
```
Read the CMakeLists.txt and debian/control in wsclean.

Tell me:
1. What are the required Build-Depends listed in CMakeLists.txt?
2. Are they all present in debian/control?
3. If anything is missing, suggest the names of the Debian packages to add.
```

**Why this works:**
The AI can read both files and cross-reference them, catching missing dependencies before the build fails.

---

### Command 3: "Full Build & Diagnostics"

**Prompt for AI:**
```
In the wsclean repo:

1. Run: debian/fetch-upstream.sh 3.6
2. Run: debuild -us -uc 2>&1 | head -100
3. If the build succeeds, run: lintian ../wsclean_3.6*.deb
4. For each warning from lintian, tell me if it's critical or can be ignored
5. Suggest the debian/* file changes needed to fix any critical warnings
```

---

## Traditional Watch File (Still Supported)

If you want to use the traditional `debian/watch` approach:

```bash
# The watch file monitors upstream:
cat debian/watch
# Output:
# version=4
# https://gitlab.com/aroffringa/wsclean/tags?sort=updated_desc .*v(\d\S+)\.tar\.gz

# Problem: uscan downloads pre-built tarballs WITHOUT submodules

# Solution: use debian/repack.sh with uscan
uscan --dehs --repack
# (The repack script adds submodules automatically)
```

---

## Git-BuildPackage Workflow (Recommended)

```bash
# Initialize the repo for gbp (one-time)
gbp clone https://salsa.debian.org/debian-astro-team/wsclean.git

# OR convert existing repo:
# (assuming you're in the wsclean repo with debian/ dir)

# Step 1: Create upstream branch with submodules
debian/fetch-upstream.sh 3.6
gbp import-orig --pristine-tar wsclean_3.6.orig.tar.gz

# Step 2: Update debian/
dch --newversion 3.6-1kern1 "New upstream v3.6"

# Step 3: Build
gbp buildpackage -us -uc

# Step 4: Push
gbp push
```

---

## File Structure After Workflow

```
/path/to/wsclean/
├── .devcontainer/
│   ├── devcontainer.json        # Codespace config
│   └── setup.sh                 # Auto-setup script
├── debian/
│   ├── fetch-upstream.sh        # Fetch + tarball creation
│   ├── repack.sh                # Enhance uscan tarballs
│   ├── update-package.sh        # Full orchestration
│   ├── changelog                # Version history
│   ├── control                  # Build depends
│   ├── watch                    # Monitor upstream (optional)
│   └── ...
├── external/
│   ├── aocommon/
│   ├── pybind11/
│   ├── radler/
│   ├── schaapcommon/
│   └── wgridder/
└── ...
```

---

## Troubleshooting

### Issue: "Submodules not included in tarball"

**Cause:** Using uscan directly without repack.sh

**Fix:**
```bash
# Use this instead:
debian/update-package.sh 3.6 noble

# OR for manual approach:
debian/fetch-upstream.sh 3.6
```

---

## For GitHub Actions / CI Integration

Add to `.github/workflows/package.yml`:

```yaml
name: Debian Package Build

on:
  push:
    tags:
      - 'v*'

jobs:
  build:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v3
        with:
          submodules: recursive
      
      - name: Install build tools
        run: sudo apt install -y build-essential debhelper devscripts cmake doxygen
      
      - name: Extract version
        id: version
        run: echo "VERSION=${GITHUB_REF##*/v}" >> $GITHUB_OUTPUT
      
      - name: Prepare upstream
        run: debian/fetch-upstream.sh ${{ steps.version.outputs.VERSION }}
      
      - name: Build package
        run: debuild -us -uc
      
      - name: Lintian check
        run: lintian ../*.deb
      
      - name: Upload artifacts
        uses: actions/upload-artifact@v3
        with:
          name: debs
          path: ../*.deb
```

---

## Key Advantages of This Workflow

✅ **Submodules included** — No more "repack" headaches  
✅ **Reproducible** — Same commands every time   
✅ **CI/CD-ready** — Can be automated in GitHub Actions, GitLab CI

---

## References

- [git-buildpackage docs](https://wiki.debian.org/UseGitBuildPackage)
- [Debian packaging guide](https://www.debian.org/doc/manuals/maint-guide/)
- [Submodule policy](https://lists.debian.org/debian-devel/2020/11/msg00180.html)
