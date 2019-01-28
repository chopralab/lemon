trap { Write-Error $_; Exit 1 }

cd c:\lemon
New-Item build -ItemType Directory
cd build

cmake .. -DLEMON_TEST_ASYNC=ON -DLEMON_LINK_SHARED=$env:SHARED -DLEMON_BUILD_PYTHON=ON -Wno-dev -G $env:generator
cmake --build . --config $env:BUILD_TYPE -- /v:q /m
