function(add_cpp_prog _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    set(_name_ lm_${_name_})

    add_executable(${_name_} ${_file_})
    set_target_properties(${_name_} PROPERTIES LINKER_LANGUAGE CXX)

    if(NOT MSVC)
        if (${COVERAGE})
            target_compile_options(${_name_} PRIVATE --coverage)
            target_link_libraries(${_name_} PRIVATE --coverage)
        endif()
    endif()

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
endfunction()

file(GLOB all_test_files
    ${CMAKE_CURRENT_SOURCE_DIR}/count/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geometry/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/interactions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/misc/*.cpp
)

foreach(test_file IN LISTS all_test_files)
    add_cpp_prog(${test_file})
endforeach(test_file)
