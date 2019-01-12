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

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
if (${LEMON_BUILD_PYTHON})
    find_package(Boost COMPONENTS system filesystem program_options python REQUIRED)
else()
    find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
endif()

include_directories(${Boost_INCLUDE_DIRS})

if (${LEMON_BUILD_DOCS})
    add_subdirectory(doc)
endif()

add_subdirectory(progs)

if (${LEMON_BUILD_PYTHON})
    find_package(PythonLibs)
    add_executable(lemon_python lang/lemon_python.cpp)

    target_include_directories(lemon_python PRIVATE
        ${PYTHON_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
    )

    target_link_libraries(lemon_python 
        ${CHEMFILES_LIBRARY}
        Boost::filesystem
        Boost::program_options
        Boost::python
        ${PYTHON_LIBRARIES}
        pthread
    )
endif()

enable_testing()
add_subdirectory(test)
