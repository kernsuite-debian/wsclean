#!/bin/bash
# WSClean Debian Packaging Workflow Status

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      WSClean AI-Augmented Debian Packaging Setup              ║"
echo "║                         ✅ READY                               ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

cd "$(dirname "$0")/.." || exit 1

echo "Package Files"
echo "├── debian/fetch-upstream.sh ........... Clone + Submodules + Tarball"
echo "├── debian/repack.sh ................... Enhance uscan Downloads"
echo "├── debian/update-package.sh ........... Full Orchestration"
echo "├── debian/inspect-deps.py ............. CMake→Debian Dependency Checker"
echo "├── debian/AI-PACKAGING.md ............. Full Documentation"
echo "├── debian/QUICKREF.md ................. Quick Reference"
echo "└── debian/watch ...................... Monitor Upstream (optional)"
echo ""

echo "Dev Container"
echo "├── .devcontainer/devcontainer.json ... Codespace Configuration"
echo "└── .devcontainer/setup.sh ............ Auto-Install Tools"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "QUICK START"
echo ""
echo "  1. Fetch upstream v3.6 with submodules:"
echo "     $ debian/fetch-upstream.sh 3.6"
echo ""
echo "  2. Or run full workflow (fetch + build check + changelog):"
echo "     $ debian/update-package.sh 3.6 noble"
echo ""
echo "  3. Build the package:"
echo "     $ debuild -us -uc"
echo ""
echo "  4. Check quality:"
echo "     $ lintian ../wsclean_3.6*.deb"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "DOCUMENTATION"
echo ""
echo "  • Full Guide:      debian/AI-PACKAGING.md"
echo "  • Quick Reference: debian/QUICKREF.md"
echo "  • For Agents:   Copy prompts from debian/QUICKREF.md"
echo ""

echo "COMMANDS"
echo ""
echo "  Prompt 1 (Update to 3.6):"
echo "    'Run debian/update-package.sh 3.6 noble, then debuild -us -uc'"
echo ""
echo "  Prompt 2 (Check Dependencies):"
echo "    'Run debian/inspect-deps.py and compare with CMakeLists.txt'"
echo ""
echo "  Prompt 3 (Full Diagnostic):"
echo "    'Build wsclean from scratch and report cmake errors'"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Check if required tools are available
echo "SYSTEM CHECK"
echo ""

missing=0
for tool in git cmake; do
    if command -v "$tool" &>/dev/null; then
        echo "   ✓ $tool"
    else
        echo "   ✗ $tool (missing, install with: sudo apt install $tool)"
        missing=$((missing+1))
    fi
done

if [ "$missing" -gt 0 ]; then
    echo ""
    echo "   Install build environment:"
    echo "   $ sudo apt install build-essential debhelper devscripts dh-make git-buildpackage"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "Next step: Run a test fetch"
echo "   $ debian/fetch-upstream.sh 3.6"
echo ""
