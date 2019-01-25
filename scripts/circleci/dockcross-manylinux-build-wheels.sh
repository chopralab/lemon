#!/bin/bash
set -e

# Pull dockcross manylinux images
docker pull dockcross/manylinux-x64

# Generate dockcross scripts
docker run dockcross/manylinux-x64 > /tmp/dockcross-manylinux-x64
chmod u+x /tmp/dockcross-manylinux-x64

pushd .
mkdir deps
cd deps
curl -L http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.gz -O
gunzip boost_1_58_0.tar.gz
tar xf boost_1_58_0.tar
cd boost_1_58_0/
sh bootstrap.sh
./bjam cxxflags=-fPIC cflags=-fPIC -a --with-filesystem --with-program_options
popd

# Build wheels
mkdir -p dist
cd lemon
DOCKER_ARGS="-v $HOME/dist:/work/dist/ -v $HOME/deps:/deps"
/tmp/dockcross-manylinux-x64 \
  -a "$DOCKER_ARGS" \
"/work/scripts/circleci/manylinux-build-wheels.sh" "$@"
