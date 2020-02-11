include(ExternalProject)

set(chemfiles_location ${CMAKE_CURRENT_BINARY_DIR}/chemfiles_install)

if (${BUILD_SHARED_LIBS})
    add_library(chemfiles SHARED IMPORTED)
    set_property(
        TARGET chemfiles
        PROPERTY IMPORTED_LOCATION
        ${chemfiles_location}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}chemfiles${CMAKE_SHARED_LIBRARY_SUFFIX}
    )

    if (MSVC)
        set_property(
            TARGET chemfiles
            PROPERTY IMPORTED_IMPLIB
            ${chemfiles_location}/lib/chemfiles.lib
        )

        set(chemfiles_byproduct ${chemfiles_location}/lib/chemfiles.lib)
    else()
        set(chemfiles_byproduct ${chemfiles_location}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}chemfiles${CMAKE_SHARED_LIBRARY_SUFFIX})
    endif()
else()
    add_library(chemfiles STATIC IMPORTED)
    set_property(
        TARGET chemfiles
        PROPERTY IMPORTED_LOCATION
            ${chemfiles_location}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}chemfiles${CMAKE_STATIC_LIBRARY_SUFFIX}
    )

    set(chemfiles_byproduct ${chemfiles_location}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}chemfiles${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

ExternalProject_Add(
    CHEMFILES
    GIT_REPOSITORY https://github.com/chemfiles/chemfiles.git
    GIT_TAG master
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/chemfiles_build
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/chemfiles_install
    CMAKE_CACHE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DCMAKE_BUILD_TYPE:STRING=Release
        -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
    EXCLUDE_FROM_ALL TRUE
    BUILD_BYPRODUCTS ${chemfiles_byproduct}
)

install(
    DIRECTORY
        ${chemfiles_location}/lib
        ${chemfiles_location}/include
    DESTINATION
        "${CMAKE_INSTALL_PREFIX}"
)

add_dependencies(chemfiles CHEMFILES)

set_target_properties(
    chemfiles PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES
        ${chemfiles_location}/include
)
