cmake_minimum_required(VERSION 3.2)

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(InstallChemfiles)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

find_package(Boost REQUIRED system filesystem)

include_directories(${Boost_INCLUDE_DIRS})

file(GLOB_RECURSE benchmarker_sources lib/**.cpp)
add_library(benchmarker STATIC ${benchmarker_sources})
add_dependencies(benchmarker chemfiles)

target_include_directories(benchmarker PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CHEMFILES_INCLUDE_DIR}
)

add_subdirectory(progs)

enable_testing()
add_subdirectory(test)
