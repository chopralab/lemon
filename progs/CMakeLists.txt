function(add_cpp_prog _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    
    add_executable(${_name_} ${_file_})
    set_target_properties(${_name_} PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(${_name_} SYSTEM PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/external
    )

    target_link_libraries( ${_name_}
        benchmarker
        chemfiles
        ${Boost_LIBRARIES}
        pthread
    )
endfunction()

file(GLOB all_test_files
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

foreach(test_file IN LISTS all_test_files)
    add_cpp_prog(${test_file})
endforeach(test_file)
