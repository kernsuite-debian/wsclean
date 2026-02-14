#!/usr/bin/env python3
"""
AI-Assistant Debian Package Inspector

This tool helps AI agents understand the wsclean build dependencies by:
1. Parsing CMakeLists.txt for find_package() calls
2. Mapping CMake modules to Debian packages
3. Comparing against debian/control
4. Suggesting missing dependencies

Usage:
  ./debian/inspect-deps.py          # Show analysis
  ./debian/inspect-deps.py --fix    # Suggest debian/control changes
"""

import re
import subprocess
from pathlib import Path

# Mapping from CMake find_package() names to Debian packages
CMAKE_TO_DEBIAN = {
    "FFTW3": "libfftw3-dev",
    "Casacore": "casacore-dev",
    "GSL": "libgsl-dev",
    "Boost": "libboost-dev",
    "CFITSIO": "libcfitsio-dev",
    "HDF5": "libhdf5-dev",
    "LAPACK": "liblapack-dev",
    "MPI": "libopenmpi-dev",
    "Python": "python3-dev",
    "Doxygen": "doxygen",
}

def parse_cmake_requirements():
    """Extract find_package() calls from CMakeLists.txt"""
    cmake_path = Path("CMakeLists.txt")
    if not cmake_path.exists():
        print("❌ CMakeLists.txt not found")
        return []
    
    content = cmake_path.read_text()
    # Match: find_package(XXX ...)
    matches = re.findall(r'find_package\s*\(\s*(\w+)', content, re.IGNORECASE)
    
    return list(set(matches))

def parse_debian_control():
    """Extract Build-Depends from debian/control"""
    control_path = Path("debian/control")
    if not control_path.exists():
        print("❌ debian/control not found")
        return []
    
    content = control_path.read_text()
    # Find Build-Depends: line and collect all package names
    match = re.search(r'Build-Depends:\s*(.*?)(?=^[A-Z]|\Z)', content, re.MULTILINE | re.DOTALL)
    
    if not match:
        return []
    
    deps_text = match.group(1)
    # Split by comma and clean up
    packages = [pkg.strip().split(' ')[0] for pkg in deps_text.split(',')]
    return [p for p in packages if p]

def main():
    print("🔍 WSClean Debian Dependency Inspector\n")
    
    cmake_deps = parse_cmake_requirements()
    debian_deps = parse_debian_control()
    
    print(f"📦 CMake find_package() calls found: {len(cmake_deps)}")
    for dep in sorted(cmake_deps):
        debian_pkg = CMAKE_TO_DEBIAN.get(dep, "?")
        status = "✓" if any(debian_pkg in d for d in debian_deps) else "✗"
        print(f"   {status} {dep:20} → {debian_pkg}")
    
    print(f"\n📚 Debian Build-Depends found: {len(debian_deps)}")
    print(f"   {', '.join(sorted(debian_deps)[:5])}...")
    
    print("\n" + "="*60)
    
    # Check for gaps
    gaps = []
    for cmake_dep in cmake_deps:
        debian_pkg = CMAKE_TO_DEBIAN.get(cmake_dep)
        if debian_pkg and not any(debian_pkg in d for d in debian_deps):
            gaps.append(debian_pkg)
    
    if gaps:
        print("\n⚠️  Potentially missing dependencies:")
        for pkg in gaps:
            print(f"   • {pkg}")
    else:
        print("\n✅ All CMake requirements appear to be in debian/control")

if __name__ == "__main__":
    main()
