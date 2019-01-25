# Find the Chemfiles library
#
# Defines:
#
#  CHEMFILES_FOUND        - system has Spglib
#  CHEMFILES_INCLUDE_DIRS - the Spglib include directories
#  CHEMFILES_LIBRARY      - The Spglib library
#
find_path(CHEMFILES_INCLUDE_DIR chemfiles.hpp HINTS ${CHEMFILES_ROOT_DIR}/include)
find_library(CHEMFILES_LIBRARY NAMES chemfiles HINTS ${CHEMFILES_ROOT_DIR}/lib)

set(CHEMFILES_INCLUDE_DIRS "${CHEMFILES_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CHEMFILES DEFAULT_MSG CHEMFILES_INCLUDE_DIR
                                  CHEMFILES_LIBRARY)

mark_as_advanced(CHEMFILES_INCLUDE_DIR CHEMFILES_LIBRARY)
