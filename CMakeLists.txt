cmake_minimum_required(VERSION 3.2)

project(lemon C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

option(LEMON_LINK_SHARED "Link lemon programs dynamically" ON)
option(LEMON_TEST_ASYNC "Build an additional test for async parallelism" OFF)
option(LEMON_BUILD_DOCS "Build documentation" OFF)
option(LEMON_BENCHMARK "Should the programs be benchmarked?" OFF)
option(LEMON_BUILD_PROGS "Should the tests and other programs be build" ON)

option(LEMON_BUILD_PYTHON "Build Python support" OFF)
option(LEMON_EXTERNAL_CHEMFILES "Should we use an external chemfiles?" OFF)

if (${LEMON_EXTERNAL_CHEMFILES})
    find_package(Chemfiles)
else()
    include(InstallChemfiles)
endif()

include(CompilerFlags)

if ((NOT ${LEMON_LINK_SHARED} OR MSVC) OR SKBUILD)
    set(Boost_USE_STATIC_LIBS ON)
endif()

set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

if (${LEMON_BUILD_DOCS})
    add_subdirectory(doc)
endif()

if (${LEMON_BUILD_PROGS})
    enable_testing()

    add_subdirectory(progs)
    add_subdirectory(test)
endif()

add_subdirectory(lang)
