if (${LEMON_BUILD_PYTHON})
    find_package(PythonLibs)
    add_executable(lemon_python ${CMAKE_CURRENT_SOURCE_DIR}/lemon_python.cpp)

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

    add_test(NAME py_simple
        COMMAND lemon_python -w ${CMAKE_SOURCE_DIR}/tests/files/rcsb_hadoop/ -p simple.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )

    add_test(NAME py_count
        COMMAND lemon_python -w ${CMAKE_SOURCE_DIR}/tests/files/rcsb_hadoop/ -p count.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )

    add_test(NAME py_invalid
        COMMAND lemon_python -w ${CMAKE_SOURCE_DIR}/tests/files/rcsb_hadoop/ -p invalid.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )

endif()
