cmake_minimum_required(VERSION 3.2)

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(InstallChemfiles)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(Boost REQUIRED system filesystem)
include_directories(${Boost_INCLUDE_DIRS})

if (${COVERAGE})
    message(STATUS "Using coverage")
endif()

enable_testing()

add_subdirectory(progs)
add_subdirectory(test)
