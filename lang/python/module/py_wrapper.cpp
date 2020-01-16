#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"
#include "lemon/geometry.hpp"
#include "lemon/tmalign.hpp"
#include "lemon/xscore.hpp"

#include "chemfiles.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-prototypes"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wmissing-braces"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated"
#pragma clang diagnostic ignored "-Wsign-conversion"

#pragma GCC diagnostic ignored "-Wold-style-cast"

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

#pragma clang diagnostic pop
