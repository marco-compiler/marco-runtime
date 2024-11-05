#!/bin/sh

src_path=$1
install_path=$2

# Create folders.
mkdir -p marco-runtime/DEBIAN
mkdir -p marco-runtime/usr/lib/marco-runtime

# Copy the control file.
cp "${src_path}/.jenkins/package/debian-12/control" marco-runtime/DEBIAN/control

# Copy the libraries.
cp "${install_path}"/lib/*.a  marco-runtime/usr/lib/marco-runtime
cp "${install_path}"/lib/*.so marco-runtime/usr/lib/marco-runtime
chmod +x marco-runtime/usr/lib/marco-runtime/*.so

# Build the package.
dpkg-deb --build marco-runtime

# Clean the work directory.
rm -rf marco-runtime

# Extract version and architecture from the control file.
version=$(grep "^Version:" "${src_path}/.jenkins/package/debian-12/control" | cut -d' ' -f2)
architecture=$(grep "^Architecture:" "${src_path}/.jenkins/package/debian-12/control" | cut -d' ' -f2)

# Rename the package.
mv marco-runtime.deb marco-runtime-${version}_${architecture}.deb
