#!/usr/bin/env bash

if [[ "${DEPLOY_MACOSX_WHEEL}" == "true" ]]; then
    exit 0
fi

cd $TRAVIS_BUILD_DIR
mkdir build && cd build

if [[ "$TRAVIS_OS_NAME" == "windows" ]]; then
    export PATH="/c/Program Files/CMake/bin":$PATH

    choco install python --version 3.6.8
    export PYTHON_EXECUTABLE="/c/Python36/python.exe"
    export PATH=$PATH:"/c/Python36/"

    # We can let Appveyor handle 32-bit for now
    export PATH=$PATH:"/c/Program Files (x86)/Microsoft Visual Studio/2017/BuildTools/MSBuild/15.0/Bin"
    export CMAKE_GENERATOR="Visual Studio 15 2017 Win64"

    export BUILD_ARGS="-verbosity:minimal -m:2"
    export BUILD_TYPE="Release"
else
    export CMAKE_GENERATOR="Unix Makefiles"
    export BUILD_ARGS="-j2"
fi

cmake .. -G "${CMAKE_GENERATOR}" \
         -DCOVERAGE=${COVERAGE} \
         -DBUILD_SHARED_LIBS=${SHARED_LIBS} \
         -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
         -DLEMON_TEST_ASYNC=${ASYNC} \
         -DLEMON_BUILD_PYTHON=ON \
         -DPYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE}

cmake --build . --config ${BUILD_TYPE} -- ${BUILD_ARGS}
ctest -j2 --output-on-failure -C ${BUILD_TYPE}
