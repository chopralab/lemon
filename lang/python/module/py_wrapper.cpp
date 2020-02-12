#define _CRT_SECURE_NO_WARNINGS

#include "lemon/lemon.hpp"
#include "lemon/external/gaurd.hpp"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wmissing-prototypes"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wmissing-braces"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated"
#pragma clang diagnostic ignored "-Wsign-conversion"
#elif __GNUC__
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif

LEMON_EXTERNAL_FILE_PUSH
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

void add_lemon_features(py::module& m);
void add_chemfiles_features(py::module& m);

// Pack the Base class wrapper into a module
PYBIND11_MODULE(lemon, m) {
    using namespace chemfiles;
    using namespace lemon;

    add_lemon_features(m);
    add_chemfiles_features(m);
}
LEMON_EXTERNAL_FILE_POP
