if (${LEMON_BUILD_PYTHON} OR SKBUILD)

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/pybind11_v2.2.4.tar.gz
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_subdirectory(
        ${CMAKE_CURRENT_BINARY_DIR}/pybind11-2.2.4/
        ${CMAKE_CURRENT_BINARY_DIR}/pybind11-2.2.4/
        EXCLUDE_FROM_ALL
    )
endif()

if (${LEMON_BUILD_PYTHON})
    add_executable(lemon_python ${CMAKE_CURRENT_SOURCE_DIR}/lemon_python.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/py_wrapper.cpp)

    target_include_directories(lemon_python PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
        ${PYTHON_INCLUDE_DIRS}
    )

    if (MSVC)
        target_link_libraries(lemon_python
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            pybind11::embed
            ws2_32
        )
    else()
        target_link_libraries(lemon_python
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            pybind11::embed
            pthread
        )
    endif()

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
endif()

if(SKBUILD)
    pybind11_add_module(lemon py_wrapper.cpp)

    target_include_directories(lemon PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
    )

    if (MSVC)
        target_link_libraries(lemon PRIVATE
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            ws2_32
        )
    else()
        target_link_libraries(lemon PRIVATE
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            pthread
        )
    endif()

    install(TARGETS lemon LIBRARY COMPONENT lemon DESTINATION "lemon_proteins")
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py COMPONENT python DESTINATION "lemon_proteins")
endif()
