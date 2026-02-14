# WSClean Debian Packaging - Quick Reference

## 🚀 One-Liner Package Update

```bash
debian/update-package.sh 3.6 noble && debuild -us -uc
```

---

## 📋 Common Tasks

### Download & Prepare Upstream (with Submodules)
```bash
debian/fetch-upstream.sh 3.6
# Output: wsclean_3.6.orig.tar.gz (in parent dir)
```

### Update Changelog
```bash
dch --newversion 3.6-1kern1 \
    --distribution noble \
    "New upstream v3.6"
```

### Build Package (unsigned, for testing)
```bash
debuild -us -uc
```

### Check Package Quality
```bash
lintian ../wsclean_3.6*.deb
```

### Check Build Dependencies
```bash
debian/inspect-deps.py
```

---

## 🤖 For AI Agents (Copy-Paste Prompts)

### Prompt 1: Upgrade to New Version
```
I need wsclean v3.6 packaged for Debian noble.

Run these commands in order and report back:
1. debian/update-package.sh 3.6 noble
2. debuild -us -uc 2>&1 | tail -20
3. lintian ../wsclean_3.6*.deb 2>&1 | head -10

For each error/warning, suggest a fix.
```

### Prompt 2: Verify Dependencies
```
Read:
- CMakeLists.txt (root level)
- debian/control

Tell me:
1. What does CMakeLists.txt require (find_package)?
2. Is each requirement in debian/control Build-Depends?
3. What packages should be added/removed?
4. Run: debian/inspect-deps.py
```

### Prompt 3: Full Diagnostic Build
```
in wsclean repo, run:

mkdir -p /tmp/wsclean-build
cd /tmp/wsclean-build
cmake /workspaces/wsclean -DCMAKE_BUILD_TYPE=Release
cmake --build . 2>&1 | grep -E 'error|warning' | head -20

If it builds, also run:
cpack --config CPackConfig.cmake

Tell me what succeeded and what failed.
```

---

## 📂 File Locations

| File | Purpose |
|------|---------|
| `debian/fetch-upstream.sh` | Clone + submodules + tarball |
| `debian/repack.sh` | Enhance uscan tarballs |
| `debian/update-package.sh` | Full orchestration (fetch + build + check) |
| `debian/inspect-deps.py` | Check build dependencies |
| `debian/AI-PACKAGING.md` | Full documentation |
| `.devcontainer/devcontainer.json` | Codespace config |
| `.devcontainer/setup.sh` | Auto-install tools |

---

## 🔧 Environment Setup

### First Time (Codespace)
```bash
# Already done! Just verify:
which debuild
which dch
which lintian
```

### Local Machine
```bash
sudo apt install build-essential debhelper devscripts dh-make \
    git-buildpackage pristine-tar cmake doxygen
```

---

## ⚠️ Troubleshooting

| Issue | Fix |
|-------|-----|
| "Submodules missing" | Use `debian/fetch-upstream.sh`, not git clone |
| "Build fails: pkg-config" | Run `debian/inspect-deps.py` to find missing -dev |
| "lintian warnings" | Read full output, not just summary |
| "Can't find CMakeLists.txt" | Must run from repo root: `cd /path/to/wsclean` |

---

## 📚 Next Steps

1. **Try it:** `debian/fetch-upstream.sh 3.7`
2. **Check it:** `ls -lh ../wsclean_3.7.orig.tar.gz`
3. **Build it:** See "Common Tasks" above
4. **Review it:** `lintian ../wsclean*.deb`

---

For detailed info, see: `debian/PACKAGING.md`
