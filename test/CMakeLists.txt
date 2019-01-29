
function(add_cpp_test _file_ _name_)    
    add_executable(${_name_} ${_file_})
    if (NOT ${LEMON_EXTERNAL_CHEMFILES})
        add_dependencies(${_name_} chemfiles)
    endif()
    set_target_properties(${_name_} PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(${_name_} SYSTEM PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/external
    )

    if(NOT MSVC)
        if (${COVERAGE})
            target_compile_options(${_name_} PRIVATE --coverage)
            target_link_libraries(${_name_} PRIVATE --coverage)
        endif()
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

    add_test(NAME ${_name_}
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${_name_}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    if (MSVC AND ${LEMON_LINK_SHARED})
        set_tests_properties( ${_name_}
            PROPERTIES ENVIRONMENT "PATH=${Boost_LIBRARY_DIRS}\;${CMAKE_BINARY_DIR}\\chemfiles\\bin"
        )
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND ${LEMON_LINK_SHARED})
        set_tests_properties( ${_name_}
            PROPERTIES ENVIRONMENT "DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/chemfiles/lib"
        )
    endif()
endfunction()

file(GLOB all_test_files
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

foreach(test_file IN LISTS all_test_files)
    get_filename_component(test_name ${test_file} NAME_WE)
    add_cpp_test(${test_file} ${test_name})
endforeach(test_file)

if(${LEMON_TEST_ASYNC})
    add_cpp_test(${CMAKE_CURRENT_SOURCE_DIR}/hadoop.cpp async)
    target_compile_definitions(async PRIVATE LEMON_USE_ASYNC=1)
endif()
