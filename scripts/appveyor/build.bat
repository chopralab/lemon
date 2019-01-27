cd C:\lemon
mkdir build
cd build
cmake .. -DLEMON_TEST_ASYNC=ON -DLEMON_LINK_SHARED=%SHARED% -DLEMON_BUILD_PYTHON=ON
cmake --build . --config %BUILD_TYPE% -- /v:q /m
