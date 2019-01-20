if (${LEMON_BUILD_PYTHON})
    find_package(PythonInterp)

    if (${PYTHONINTERP_FOUND})
        if (UNIX AND NOT APPLE)
            if (${PYTHON_VERSION_MAJOR} EQUAL 3)
                find_package(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX})
                find_package(PythonLibs 3 REQUIRED)
            else()
                find_package(Boost COMPONENTS python)
                find_package(PythonLibs REQUIRED)
            endif()
        else()
            if (${PYTHON_VERSION_MAJOR} EQUAL 3)
                find_package(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
                find_package(PythonLibs 3 REQUIRED)
            else()
                find_package(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
                find_package(PythonLibs REQUIRED)
            endif()
        endif()
    else()
        message("Python not found")
    endif()

    python_add_module(lemon py_wrapper.cpp)
    target_include_directories(lemon PRIVATE
        ${PYTHON_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
    )

    add_executable(lemon_python ${CMAKE_CURRENT_SOURCE_DIR}/lemon_python.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/py_wrapper.cpp)

    target_include_directories(lemon_python PRIVATE
        ${PYTHON_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include
        ${CHEMFILES_INCLUDE_DIR}
    )

    if (UNIX AND NOT APPLE)
        target_link_libraries(lemon
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            Boost::python
            ${PYTHON_LIBRARIES}
            pthread
        )

        target_link_libraries(lemon_python
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            Boost::python
            ${PYTHON_LIBRARIES}
            pthread
        )
    elseif(WIN32)
        target_link_libraries(lemon
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            Boost::python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}
            ${PYTHON_LIBRARIES}
        )

        target_link_libraries(lemon_python
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            Boost::python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}
            ${PYTHON_LIBRARIES}
        )
    else()
        target_link_libraries(lemon
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            Boost::python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}
            ${PYTHON_LIBRARIES}
            pthread
        )

        target_link_libraries(lemon_python
            ${CHEMFILES_LIBRARY}
            Boost::filesystem
            Boost::program_options
            Boost::python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}
            ${PYTHON_LIBRARIES}
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
