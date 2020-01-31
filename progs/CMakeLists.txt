function(setup_prog _name_)
    if(${LEMON_TEST_ASYNC})
        target_compile_definitions(${_name_} PRIVATE LEMON_USE_ASYNC=1)
    endif()

    if (${LEMON_BENCHMARK})
        target_compile_definitions(${_name_} PRIVATE LEMON_BENCHMARK=1)
    endif()

    target_link_libraries( ${_name_} PUBLIC lemon)
    if(NOT MSVC)
        if (APPLE)
            set_target_properties(${_name_} PROPERTIES
                MACOSX_RPATH 1
                INSTALL_RPATH "@loader_path/../lib"
            )
        else()
            set_target_properties(${_name_} PROPERTIES
                INSTALL_RPATH "\$ORIGIN/../lib"
            )
        endif()
    endif()

    install( TARGETS ${_name_} RUNTIME DESTINATION bin)
endfunction(setup_prog)

function(add_cpp_prog _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    set(_name_ lm_${_name_})

    add_executable(${_name_} ${_file_})
    set_target_properties(${_name_} PROPERTIES LINKER_LANGUAGE CXX)

    setup_prog(${_name_})
endfunction(add_cpp_prog)

file(GLOB all_prog_files
    ${CMAKE_CURRENT_SOURCE_DIR}/count/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geometry/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/interactions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/misc/*.cpp
)

foreach(prog_file IN LISTS all_prog_files)
    add_cpp_prog(${prog_file})
endforeach(prog_file)

# Add DUBS

add_executable( dubs dubs/dubs.cpp dubs/parser.cpp )

target_include_directories(dubs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/dubs)
setup_prog(dubs)
