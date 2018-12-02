cmake_minimum_required(VERSION 3.2)

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(InstallChemfiles)

option(LEMON_TEST_ASYNC "Build an additional test for async parallelism" OFF)
option(LEMON_BUILD_DOCS "Build documentation" OFF)

if (${LEMON_TEST_ASYNC})
    set(CMAKE_CXX_STANDARD 14)
else()
    set(CMAKE_CXX_STANDARD 11)
endif()

set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "Compiler")
set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER} CACHE FILEPATH "Compiler")

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

if (${COVERAGE})
    message(STATUS "Using coverage")
endif()

if (${LEMON_BUILD_DOCS})
    add_subdirectory(doc)
endif()

add_subdirectory(progs)

enable_testing()
add_subdirectory(test)
