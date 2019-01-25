#!/usr/bin/env bash
set -e 

SCRIPT_DIR=$(cd $(dirname $0) || exit 1; pwd)

# Install python versions
source "${SCRIPT_DIR}/macpython-install-python.sh"

MACPYTHON_PY_PREFIX=/Library/Frameworks/Python.framework/Versions

# Versions can be restricted by passing them in as arguments to the script
# For example,
# macpython-build-wheels.sh 2.7 3.5
if [[ $# -eq 0 ]]; then
  PYBINARIES=(${MACPYTHON_PY_PREFIX}/*)
else
  PYBINARIES=()
  for version in "$@"; do
    PYBINARIES+=(${MACPYTHON_PY_PREFIX}/*${version}*)
  done
fi

VENVS=()
mkdir -p ${SCRIPT_DIR}/../venvs
for PYBIN in "${PYBINARIES[@]}"; do
    if [[ $(basename $PYBIN) = "Current" ]]; then
      continue
    fi
    py_mm=$(basename ${PYBIN})
    VENV=${SCRIPT_DIR}/../venvs/${py_mm}
    VENVS+=(${VENV})
done

# Since the python interpreter exports its symbol (see [1]), python
# modules should not link against any python libraries.
# To ensure it is not the case, we configure the project using an empty
# file as python library.
#
# [1] "Note that libpythonX.Y.so.1 is not on the list of libraries that
# a manylinux1 extension is allowed to link to. Explicitly linking to
# libpythonX.Y.so.1 is unnecessary in almost all cases: the way ELF linking
# works, extension modules that are loaded into the interpreter automatically
# get access to all of the interpreter's symbols, regardless of whether or
# not the extension itself is explicitly linked against libpython. [...]"
#
# Source: https://www.python.org/dev/peps/pep-0513/#libpythonx-y-so-1
PYTHON_LIBRARY=$(cd $(dirname $0); pwd)/libpython-not-needed-symbols-exported-by-interpreter
touch ${PYTHON_LIBRARY}

# -----------------------------------------------------------------------
# Remove previous virtualenv's
rm -rf ${SCRIPT_DIR}/../venvs
# Create virtualenv's
VENVS=()
mkdir -p ${SCRIPT_DIR}/../venvs
for PYBIN in "${PYBINARIES[@]}"; do
    if [[ $(basename $PYBIN) = "Current" ]]; then
      continue
    fi
    py_mm=$(basename ${PYBIN})
    VENV=${SCRIPT_DIR}/../venvs/${py_mm}
    VIRTUALENV_EXECUTABLE=${PYBIN}/bin/virtualenv
    ${VIRTUALENV_EXECUTABLE} ${VENV}
    VENVS+=(${VENV})
done

VENV="${VENVS[0]}"
PYTHON_EXECUTABLE=${VENV}/bin/python
${PYTHON_EXECUTABLE} -m pip install --no-cache cmake
CMAKE_EXECUTABLE=${VENV}/bin/cmake
${PYTHON_EXECUTABLE} -m pip install --no-cache ninja
NINJA_EXECUTABLE=${VENV}/bin/ninja
${PYTHON_EXECUTABLE} -m pip install --no-cache delocate
DELOCATE_LISTDEPS=${VENV}/bin/delocate-listdeps
DELOCATE_WHEEL=${VENV}/bin/delocate-wheel

# Include deps
pushd .
cd $TRAVIS_BUILD_DIR/../
mkdir deps
cd deps
curl -L http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.gz -O
gunzip boost_1_58_0.tar.gz
tar xf boost_1_58_0.tar
cd boost_1_58_0/
sh bootstrap.sh
./bjam cxxflags=-fPIC cflags=-fPIC -a --with-filesystem --with-program_options
popd

for VENV in "${VENVS[@]}"; do
    py_mm=$(basename ${VENV})
    export PYTHON_EXECUTABLE=${VENV}/bin/python
    PYTHON_INCLUDE_DIR=$( find -L ${MACPYTHON_PY_PREFIX}/${py_mm}/include -name Python.h -exec dirname {} \; )

    echo ""
    echo "PYTHON_EXECUTABLE:${PYTHON_EXECUTABLE}"
    echo "PYTHON_INCLUDE_DIR:${PYTHON_INCLUDE_DIR}"
    echo "PYTHON_LIBRARY:${PYTHON_LIBRARY}"

    # Install dependencies
    ${PYTHON_EXECUTABLE} -m pip install --upgrade -r ${SCRIPT_DIR}/../../requirements-dev.txt

    build_type="MinSizeRel"
    plat_name="macosx-10.9-x86_64"
    osx_target="10.9"

    # Clean up previous invocations
    rm -rf _skbuild
    
    # Generate wheel
    ${PYTHON_EXECUTABLE} setup.py bdist_wheel --build-type ${build_type} --plat-name ${plat_name} -G 'Unix Makefiles' -- \
      -DBOOST_ROOT:PATH=$TRAVIS_BUILD_DIR/../deps/boost_1_58_0/
      -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=${osx_target} \
      -DCMAKE_OSX_ARCHITECTURES:STRING=x86_64 \
      -DPYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE} \
      -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR} \
      -DPYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY} 
    # Cleanup
    ${PYTHON_EXECUTABLE} setup.py clean
done

#${DELOCATE_LISTDEPS} ${SCRIPT_DIR}/../dist/*.whl # lists library dependencies
#${DELOCATE_WHEEL} ${SCRIPT_DIR}/../dist/*.whl # copies library dependencies into wheel
