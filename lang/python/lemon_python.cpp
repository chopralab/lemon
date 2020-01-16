#include <string>
#include <iostream>
#include <mutex>

#pragma clang diagnostic ignored "-Wmissing-prototypes"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wmissing-braces"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated"
#pragma clang diagnostic ignored "-Wsign-conversion"

#pragma GCC diagnostic ignored "-Wold-style-cast"

// Mac OSX problems with a tolower macro
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"
#include <pybind11/embed.h>

namespace python = pybind11;

#if PY_MAJOR_VERSION >= 3
#   define INIT_MODULE PyInit_lemon
    extern "C" PyObject* INIT_MODULE();
#else
#   define INIT_MODULE initlemon
    extern "C" void INIT_MODULE();
#endif

int main(int argc, char *argv[]) {
    lemon::Options o;
    std::string py_script("lemon.py");
    o.add_option("--py_script,-p", py_script, "Python script to load");
    o.parse_command_line(argc, argv);

    // Register the module with the interpreter
    if (PyImport_AppendInittab("lemon", INIT_MODULE) == -1) {
        std::cerr << "Failed to embed lemon in to builtin modules" << std::endl;
        return 1;
    }

    python::scoped_interpreter guard{};

    auto locals = python::dict();

    try {
        python::exec("LEMON_HADOOP_DIR='" + o.work_dir() + "'\n");
        python::exec("LEMON_NUM_THREADS=" + std::to_string(o.ncpu()));
        python::eval_file(py_script, python::globals(), locals);
    } catch (python::error_already_set& err) {
        std::cerr << err.what() << std::endl;
        return 1;
    } catch (std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }

    return 0;
}
