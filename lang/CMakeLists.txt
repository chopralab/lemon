if (${LEMON_BUILD_PYTHON})

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/pybind11_v2.2.4.tar.gz
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_subdirectory(
        ${CMAKE_CURRENT_BINARY_DIR}/pybind11-2.2.4/
        ${CMAKE_CURRENT_BINARY_DIR}/pybind11-2.2.4/
        EXCLUDE_FROM_ALL
    )

    #python_add_module(_lemon py_wrapper.cpp)
    #target_include_directories(_lemon PRIVATE
    #    ${PYTHON_INCLUDE_DIRS}
    #    ${CMAKE_SOURCE_DIR}/include
    #    ${CHEMFILES_INCLUDE_DIR}
    #)

    add_executable(lemon_python ${CMAKE_CURRENT_SOURCE_DIR}/lemon_python.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/py_wrapper.cpp)

    target_include_directories(lemon_python PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
        ${PYTHON_INCLUDE_DIRS}
    )

    target_link_libraries(lemon_python
        ${CHEMFILES_LIBRARY}
        Boost::filesystem
        Boost::program_options
        pybind11::embed
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

    add_test(NAME py_small_heme
        COMMAND lemon_python -w ${CMAKE_SOURCE_DIR}/tests/files/rcsb_hadoop/ -p small_heme.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )

    add_test(NAME py_tmscore
        COMMAND lemon_python -w ${CMAKE_SOURCE_DIR}/tests/files/rcsb_hadoop/ -p tmscore.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )

    if (MSVC AND ${LEMON_LINK_SHARED})
        set_tests_properties( py_simple py_count py_invalid py_small_heme py_tmscore
            PROPERTIES ENVIRONMENT "PATH=${PATH_STRING}\;${Boost_LIBRARY_DIRS}\;${CMAKE_BINARY_DIR}\\chemfiles\\bin"
        )
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND ${LEMON_LINK_SHARED})
        set_tests_properties( py_simple py_count py_invalid py_small_heme py_tmscore
            PROPERTIES ENVIRONMENT "DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/chemfiles/lib"
        )
    endif()

    if(SKBUILD)
        message(${CMAKE_INSTALL_PREFIX})
        install(TARGETS _lemon LIBRARY DESTINATION lemon RENAME lemon)
    endif()
endif()
