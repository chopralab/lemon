#!/usr/bin/env bash

# To run it, you will need to generate a compilation database first. CMake can
# generate one for you using `cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON`. It is
# better to use the same version of clang as clang-tidy when configuring with
# cmake.
#
# CMake will generate a `compile_commands.json` file, that you should copy/link
# in the root of chemfiles sources. You can then run this script from the root
# as ./scripts/tidy.sh

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../" && pwd)

if [[ ! -f "$ROOT/compile_commands.json" ]]; then
    echo "Missing compile_commands.json. See the script for how to generate it"
    exit 1
fi

HEADERS='include/lemon/*.hpp'

pushd $ROOT > /dev/null

for file in $(find progs -name '*.?pp')
do
    echo "===== checking ${file} ======"
    clang-tidy-9 -p . $file -header-filter='.*'
done

popd > /dev/null

