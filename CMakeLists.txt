cmake_minimum_required(VERSION 3.2)

project(lemon C CXX)

set(LEMON_VERSION 0.1.2)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

option(BUILD_SHARED_LIBS "Link lemon programs dynamically" ON)
option(LEMON_TEST_ASYNC "Build an additional test for async parallelism" OFF)
option(LEMON_BUILD_DOCS "Build documentation" OFF)
option(LEMON_BENCHMARK "Should the programs be benchmarked?" OFF)
option(LEMON_BUILD_PROGS "Should the tests and other programs be build" ON)

option(LEMON_BUILD_PYTHON "Build Python support" OFF)
option(LEMON_EXTERNAL_CHEMFILES "Should we use an external chemfiles?" OFF)

if (${LEMON_EXTERNAL_CHEMFILES})
    find_package(chemfiles REQUIRED)
else()
    include(InstallChemfiles)
endif()

include(CompilerFlags)

#----------------------------------------------------------------------------------------#
# Installation configuration
#----------------------------------------------------------------------------------------#

add_library(lemon INTERFACE)
install(TARGETS lemon EXPORT lemon-targets
    INCLUDES DESTINATION include
)

target_include_directories(lemon INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

if (MSVC)
    target_link_libraries(lemon INTERFACE chemfiles ws2_32)
else()
    target_link_libraries(lemon INTERFACE chemfiles pthread)
endif()

install(DIRECTORY include/ DESTINATION include)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/lemon-config-version.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/lemon-config-version.cmake"
    @ONLY
)

# Install package files
set(PACKAGE_LOCATION lib${LIB_SUFFIX}/cmake/lemon)
install(EXPORT lemon-targets
  FILE lemon-targets.cmake
  DESTINATION ${PACKAGE_LOCATION}
)
install(FILES
    cmake/lemon-config.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/lemon-config-version.cmake"
    DESTINATION ${PACKAGE_LOCATION}
    COMPONENT Devel
)

# uninstall target
configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/uninstall.in.cmake"
    "${PROJECT_BINARY_DIR}/uninstall.cmake"
    @ONLY
)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BINARY_DIR}/uninstall.cmake
)

if (${LEMON_BUILD_DOCS})
    add_subdirectory(doc)
endif()

if (${LEMON_BUILD_PROGS})
    enable_testing()

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/launch_lemon.pbs.in"
        "${CMAKE_CURRENT_BINARY_DIR}/launch_lemon.pbs"
        @ONLY
    )

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/launch_lemon.pbs" DESTINATION bin)

    add_subdirectory(progs)
    add_subdirectory(test)
endif()

add_subdirectory(lang)

