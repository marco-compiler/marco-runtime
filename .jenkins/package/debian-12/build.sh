#!/bin/sh

# The path where the project was cloned.
src_path=$1

# The path where the project was installed.
install_path=$2

# Extract version and architecture from the control file.
version=$(grep "^Version:" "${src_path}/.jenkins/package/debian-12/control" | cut -d' ' -f2)
architecture=$(grep "^Architecture:" "${src_path}/.jenkins/package/debian-12/control" | cut -d' ' -f2)

# Create folders.
package_name=marco-runtime-${version}_${architecture}

mkdir -p "${package_name}"/DEBIAN
mkdir -p "${package_name}"/usr/lib/marco-runtime

# Copy the control file.
cp "${src_path}/.jenkins/package/debian-12/control" "${package_name}"/DEBIAN/control

# Copy the libraries.
cp "${install_path}"/lib/*.a "${package_name}"/usr/lib/marco-runtime
cp "${install_path}"/lib/*.so "${package_name}"/usr/lib/marco-runtime
chmod +x "${package_name}"/usr/lib/marco-runtime/*.so

# Build the package.
dpkg-deb --build "${package_name}"

# Clean the work directory.
rm -rf "${package_name}"
