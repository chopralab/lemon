#!/usr/bin/env bash

echo $PYPI_PASSWORD

if [[ "${DEPLOY_MACOSX_WHEEL}" == "true" ]]; then
    exit 0
fi

${CC} --version
${CXX} --version
cd $TRAVIS_BUILD_DIR
mkdir build && cd build

cmake .. -DCOVERAGE=${COVERAGE} \
         -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
         -DLEMON_TEST_ASYNC=${ASYNC} \
         -DLEMON_BUILD_PYTHON=ON \
         -DPYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE}

make -j2
ctest -j2 --output-on-failure
