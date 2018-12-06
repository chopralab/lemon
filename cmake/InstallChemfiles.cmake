include(ExternalProject)

ExternalProject_Add( CHEMFILES
    GIT_REPOSITORY https://github.com/frodofine/chemfiles.git
    GIT_TAG read_from_memory
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/chemfiles_build
    CMAKE_CACHE_ARGS    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/chemfiles
                        -DBUILD_SHARED_LIBS:BOOL=${LEMON_LINK_SHARED}
                        -DCMAKE_BUILD_TYPE:STRING=Release
)

if (${LEMON_LINK_SHARED})
    add_library(chemfiles SHARED IMPORTED)
    set_property(TARGET chemfiles PROPERTY IMPORTED_LOCATION
        ${CMAKE_CURRENT_BINARY_DIR}/chemfiles/lib/${CMAKE_SHARED_LIBRARY_PREFIX}chemfiles${CMAKE_SHARED_LIBRARY_SUFFIX}
    )
    set_property(TARGET chemfiles PROPERTY IMPORTED_IMPLIB
        ${CMAKE_CURRENT_BINARY_DIR}/chemfiles/lib/chemfiles.lib
    )
else()
    add_library(chemfiles STATIC IMPORTED)
    set_property(TARGET chemfiles PROPERTY IMPORTED_LOCATION
        ${CMAKE_CURRENT_BINARY_DIR}/chemfiles/lib/${CMAKE_STATIC_LIBRARY_PREFIX}chemfiles${CMAKE_STATIC_LIBRARY_SUFFIX}
    )
endif()

add_dependencies(chemfiles CHEMFILES)

set(CHEMFILES_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/chemfiles/include)
set(CHEMFILES_LIBRARY chemfiles)
