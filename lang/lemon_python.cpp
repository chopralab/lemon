#include <string>
#include <iostream>

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
    std::string py_derive("MyWorkflow");
    o.add_option("--py_script,-p", py_script, "Python script to load");
    o.add_option("--py_class,-c", py_derive, "Class deriving from Workflow");
    o.parse_command_line(argc, argv);

    // Register the module with the interpreter
    if (PyImport_AppendInittab("lemon", INIT_MODULE) == -1) {
        std::cerr << "Failed to embed lemon in to builtin modules" << std::endl;
        return 1;
    }

    python::scoped_interpreter guard{};

    auto locals = python::dict();
    try {
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

    // Obtained derived class from python
    python::object PythonDerived = locals[py_derive.c_str()];
    python::object py_base = PythonDerived();
    lemon::LemonPythonBase& py = py_base.cast<lemon::LemonPythonBase&>();

    auto worker = [&py](chemfiles::Frame entry, const std::string& pdbid) {
        try {
            return py.worker(&entry, pdbid);
        } catch (python::error_already_set& err) {
            std::cerr << pdbid + " " + err.what() + "\n";
        } catch (python::cast_error& err) {
            std::cerr << pdbid + " Problem with type: " + err.what() + "\n";
        } catch (std::exception& err) {
            std::cerr << pdbid + " " + err.what() + "\n";
        } catch (...) {
            std::cerr << pdbid + " unknown error." + "\n";
        }
        std::terminate();
    };

    auto collector = lemon::print_combine(std::cout);
    lemon::launch(o, worker, collector);

    try {
        py.finalize();
    } catch(python::error_already_set& err) {
        std::cerr << err.what() << std::endl;
    }

    return 0;
}
