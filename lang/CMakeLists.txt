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

function(add_python_test _name_)
    add_test(NAME py_${_name_}
        COMMAND lemon_python -w ${CMAKE_SOURCE_DIR}/test/files/rcsb_hadoop/ -p ${_name_}.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )

    if (MSVC AND ${BUILD_SHARED_LIBS})
        set_tests_properties( py_${_name_}
            PROPERTIES ENVIRONMENT "PATH=${PATH_STRING}\;${CMAKE_BINARY_DIR}/chemfiles/bin/"
        )
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND ${BUILD_SHARED_LIBS})
        set_tests_properties( py_${_name_}
            PROPERTIES ENVIRONMENT "DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/chemfiles/lib"
        )
    endif()
endfunction()

if (${LEMON_BUILD_PYTHON})
    add_executable(lemon_python ${CMAKE_CURRENT_SOURCE_DIR}/lemon_python.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/py_wrapper.cpp)

    target_include_directories(lemon_python PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${PYTHON_INCLUDE_DIRS}
    )

    target_link_libraries(lemon_python lemon pybind11::embed)

    add_python_test(simple)
    add_python_test(count)
    add_python_test(metal_ions)
    add_python_test(residues)
    add_python_test(small_molecule)
    add_python_test(small_heme)
    add_python_test(dry_small_molecule)
    add_python_test(protein_dihedral)
    add_python_test(tmscore)
    add_python_test(vina)
endif()

if(SKBUILD)
    pybind11_add_module(pylemon py_wrapper.cpp)

    target_include_directories(pylemon PRIVATE
        ${CMAKE_SOURCE_DIR}/include
    )

    target_link_libraries(pylemon lemon)

    install(TARGETS lemon LIBRARY COMPONENT lemon DESTINATION "candiy_lemon" RENAME lemon)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py COMPONENT python DESTINATION "candiy_lemon")
endif()
