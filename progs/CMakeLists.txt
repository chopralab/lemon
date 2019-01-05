function(add_cpp_prog _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    
    add_executable(${_name_} ${_file_})
    add_dependencies(${_name_} chemfiles)
    set_target_properties(${_name_} PROPERTIES LINKER_LANGUAGE CXX)

    if(MSVC)
        target_compile_options(${_name_} PRIVATE /W4)
    else()
        target_compile_options(${_name_} PRIVATE -Wall -Wextra -pedantic)
        if (${COVERAGE})
            target_compile_options(${_name_} PRIVATE --coverage)
            target_link_libraries(${_name_} PRIVATE --coverage)
        endif()
    endif()

    target_include_directories(${_name_} SYSTEM PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/external
    )

    if(${LEMON_TEST_ASYNC})
        target_compile_definitions(${_name_} PRIVATE LEMON_USE_ASYNC=1)
    endif()

    if (MSVC)
        target_link_libraries( ${_name_} PUBLIC
            ${CHEMFILES_LIBRARY}
            ws2_32
            Boost::filesystem
            Boost::program_options
        )
    else()
        target_link_libraries( ${_name_} PUBLIC
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            pthread
        )    
    endif()
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
