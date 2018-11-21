cmake_minimum_required(VERSION 3.2)

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(InstallChemfiles)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "Compiler")
set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER} CACHE FILEPATH "Compiler")

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

if (${COVERAGE})
    message(STATUS "Using coverage")
endif()

add_subdirectory(progs)

enable_testing()
add_subdirectory(test)
