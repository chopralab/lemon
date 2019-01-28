# Taken from Chemfiles

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (${LEMON_TEST_ASYNC})
    set(CMAKE_CXX_STANDARD 14)
else()
    set(CMAKE_CXX_STANDARD 11)
endif()

set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "Compiler")
set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER} CACHE FILEPATH "Compiler")

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)

set(CMAKE_REQUIRED_QUIET YES)

if(MSVC)
    add_definitions("/D NOMINMAX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
    set(CMAKE_SHARED_LINKER_FLAGS "/SUBSYSTEM:CONSOLE")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

macro(add_warning_flag _flag_)
    CHECK_CXX_COMPILER_FLAG("${_flag_}" CXX_SUPPORTS${_flag_})
    CHECK_C_COMPILER_FLAG("${_flag_}" CC_SUPPORTS${_flag_})
    if(CXX_SUPPORTS${_flag_})
        set(LEMON_CXX_WARNINGS "${LEMON_CXX_WARNINGS} ${_flag_}")
    endif()

    if(CC_SUPPORTS${_flag_})
        set(LEMON_C_WARNINGS "${LEMON_C_WARNINGS} ${_flag_}")
    endif()
endmacro()

set(LEMON_CXX_WARNINGS "")
set(LEMON_C_WARNINGS "")

macro(remove_msvc_warning _warn_)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd${_warn_}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd${_warn_}")
endmacro()

macro(remove_intel_warning _warn_)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable ${_warn_}")
endmacro()

if(MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/Wall" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
    endif()

    if(CMAKE_C_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/Wall" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Wall")
    endif()

    # Disable other warnings
    remove_msvc_warning(4061) # enumerator in switch of enum is not explicitly handled by a case label
    remove_msvc_warning(4275) # non-dll-export interface as base class of dll-export class
    remove_msvc_warning(4251) # class <> needs to be dll-export
    remove_msvc_warning(4514) # unreferenced inline function has been removed
    remove_msvc_warning(4571) # SEH exceptions no longer caught
    remove_msvc_warning(4582) # constructor is not implicitly called
    remove_msvc_warning(4583) # destructor is not implicitly called
    remove_msvc_warning(4623) # default constructor was implicitly defined as deleted
    remove_msvc_warning(4625) # copy constructor was implicitly defined as deleted
    remove_msvc_warning(4626) # assignment operator was implicitly defined as deleted
    remove_msvc_warning(4668) # not defined preprocessor macro, replacing with '0' for '#if/#elif'
    remove_msvc_warning(4627) # move assignment operator was implicitly defined as deleted
    remove_msvc_warning(4710) # function not inlined
    remove_msvc_warning(4711) # function selected for automatic inlining
    remove_msvc_warning(4774) # arguments to print_s not string literals
    remove_msvc_warning(4820) # padding added
    remove_msvc_warning(5026) # move constructor was implicitly defined as deleted
    remove_msvc_warning(5027) # move assignment operator was implicitly defined as deleted
else()
    # Add some warnings in debug mode
    # Basic set of warnings
    add_warning_flag("-Wall")
    add_warning_flag("-Wextra")
    # Initialization and convertion values
    add_warning_flag("-Wuninitialized")
    add_warning_flag("-Wconversion")
    add_warning_flag("-Wsign-conversion")
    add_warning_flag("-Wsign-promo")
    # C++11 functionalities
    add_warning_flag("-Wsuggest-override")
    add_warning_flag("-Wsuggest-final-types")
    # C++ standard conformance
    add_warning_flag("-Wpedantic")
    add_warning_flag("-pedantic")
    # The compiler is your friend
    add_warning_flag("-Wdocumentation")
    add_warning_flag("-Wdeprecated")
    add_warning_flag("-Wextra-semi")
    add_warning_flag("-Wnon-virtual-dtor")
    add_warning_flag("-Wold-style-cast")
    add_warning_flag("-Wcast-align")
    add_warning_flag("-Wunused")
    add_warning_flag("-Woverloaded-virtual")
    add_warning_flag("-Wundefined-func-template")
    add_warning_flag("-Wmissing-prototypes")
    add_warning_flag("-Wmissing-variable-declarations")

    # Disable some strong warning with clang that are OK here
    add_warning_flag("-Wno-unknown-pragmas")
    add_warning_flag("-Wno-weak-vtables")
    add_warning_flag("-Wno-weak-template-vtables")
    add_warning_flag("-Wno-switch-enum")
    # We are not doing C here
    add_warning_flag("-Wno-padded")
    # Sometime this OK
    add_warning_flag("-Wno-float-equal")
    add_warning_flag("-Wno-double-promotion")
    # Yes, chemfiles uses globals
    add_warning_flag("-Wno-exit-time-destructors")
    add_warning_flag("-Wno-global-constructors")
    # Not everyone is as smart as clang for code reachability
    add_warning_flag("-Wno-covered-switch-default")
    add_warning_flag("-Wno-unreachable-code-break")

    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        # external function definition with no prior declaration
        remove_intel_warning(1418)
        # Intel compiler is too strict in errors about 'explicit' keyword
        remove_intel_warning(2304)
        remove_intel_warning(2305)
        # parameter "args" was never referenced in variadic templates
        remove_intel_warning(869)
        # exception specification for implicitly declared virtual function ...
        # This is an issue with compiler generated destructors and noexcept
        remove_intel_warning(811)
        remove_intel_warning(809)
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LEMON_CXX_WARNINGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LEMON_C_WARNINGS}")
