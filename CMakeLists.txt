cmake_minimum_required(VERSION 3.2)

project(lemon C CXX)

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

option(LEMON_LINK_SHARED "Link lemon programs dynamically" ON)
option(LEMON_TEST_ASYNC "Build an additional test for async parallelism" OFF)
option(LEMON_BUILD_DOCS "Build documentation" OFF)

option(LEMON_BUILD_PYTHON "Build Python support" OFF)

include(InstallChemfiles)
include(CompilerFlags)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

if (${LEMON_BUILD_DOCS})
    add_subdirectory(doc)
endif()

enable_testing()

add_subdirectory(progs)
add_subdirectory(test)

add_subdirectory(lang)
