#include <string>
#include <iostream>

// Mac OSX problems with a tolower macro
#include "lemon/lemon.hpp"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

namespace python = boost::python;

#if PY_MAJOR_VERSION >= 3
#   define INIT_MODULE PyInit_lemon
    extern "C" PyObject* INIT_MODULE();
#else
#   define INIT_MODULE initlemon
    extern "C" void INIT_MODULE();
#endif

char* get_python_error() {
    PyObject *ptype, *pvalue, *ptraceback;
    PyErr_Fetch(&ptype, &pvalue, &ptraceback);
    #if PY_MAJOR_VERSION >= 3
        PyObject* pyStr = PyUnicode_AsEncodedString(pvalue, "utf-8", "Error ~");
        return PyBytes_AS_STRING(pyStr);
    #else
        return PyString_AsString(pvalue);
    #endif
}

int main(int argc, char *argv[]) {
    lemon::Options o;
    std::string py_script("lemon.py");
    std::string py_derive("MyWorkflow");
    o.add_option("py_script,p", py_script, "Python script to load");
    o.add_option("py_class,c", py_script, "Class deriving from Workflow");
    o.parse_command_line(argc, argv);

    Py_Initialize();

    // Register the module with the interpreter
    if (PyImport_AppendInittab("lemon", INIT_MODULE) == -1) {
        std::cerr << "Failed to embed lemon in to builtin modules" << std::endl;
        return 1;
    }

    // Retrieve the main module
    python::object main = python::import("__main__");
  
    // Retrieve the main module's namespace
    python::object global(main.attr("__dict__"));

    // Get around a Boost bug in exec_file
    std::ifstream t(py_script);
    std::string str((std::istreambuf_iterator<char>(t)),
                     std::istreambuf_iterator<char>());
    try {
        python::exec(str.c_str(), global, global);
    } catch (python::error_already_set const &) {
        std::cerr << "Error in '" << py_script <<
                     "': " << get_python_error() << "\n";
        return 1;
    }

    // Obtained derived class from python
    python::object PythonDerived = global[py_derive.c_str()];
    python::object py_base = PythonDerived();
    lemon::LemonPythonBase& py = python::extract<lemon::LemonPythonBase&>(py_base);

    auto worker = [&py](chemfiles::Frame complex, const std::string& pdbid) {
        try {
            return py.worker(complex, pdbid);
        } catch (...) {
            return pdbid + " " + get_python_error() + "\n";
        }
    };

    lemon::launch<lemon::print_combine>(o, worker, std::cout);

    try {
        py.finalize();
    } catch(...) {
        // Ignore, the user probably didn't define the finalize function.
    }

    return 0;
}
