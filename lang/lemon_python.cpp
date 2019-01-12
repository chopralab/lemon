#include <string>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "lemon/lemon.hpp"

namespace python = boost::python;

class LemonPythonBase : public boost::noncopyable {
public:
    virtual ~LemonPythonBase() {};
    virtual std::string worker(chemfiles::Frame& frame) = 0;
};

struct LemonPythonWrap : LemonPythonBase, python::wrapper<LemonPythonBase> {
    virtual std::string worker(chemfiles::Frame& frame) override {
        return this->get_override("worker")(boost::ref(frame));
    }
};

template<typename T, typename ret, int maxval>
ret get_index(const T& v, int i) {
    if (i >= 0) {
        if (i >= maxval) {
            PyErr_SetString(PyExc_IndexError, "index is too large");
            python::throw_error_already_set();
        }
        return v[static_cast<size_t>(i)];
    }
    if (-i > maxval) {
        PyErr_SetString(PyExc_IndexError, "index is too small");
        python::throw_error_already_set();
    }
    return v[static_cast<size_t>(maxval + i)];
}

template<typename T>
bool check(const chemfiles::optional<T>& o) {
    return o != chemfiles::nullopt;
}

template<typename T>
T get(const chemfiles::optional<T>& o) {
    if (o) {
        return *o;
    } else {
        PyErr_SetString(PyExc_IndexError, "Cannot dereference nullopt");
        python::throw_error_already_set();
    }
    chemfiles::unreachable();
}

template<typename T>
chemfiles::optional<const chemfiles::Property&> getp(const T& container,
                                                     std::string name) {
    return container.get(name);
}

// Pack the Base class wrapper into a module
BOOST_PYTHON_MODULE(lemon) {
    using namespace chemfiles;
    using boost::noncopyable;

    python::class_<LemonPythonWrap, noncopyable>("Workflow");

    /**************************************************************************
     * Optional
     **************************************************************************/
    python::class_<optional<uint64_t>> ("OptionalUInt64", python::no_init)
        .def("check", check<uint64_t>)
        .def("get", get<uint64_t>);

    python::class_<optional<std::string>> ("OptionalString", python::no_init)
        .def("check", check<std::string>)
        .def("get", get<std::string>);

    python::class_<optional<double>> ("OptionalDouble", python::no_init)
        .def("check", check<double>)
        .def("get", get<double>);

    python::class_<optional<const Property&>> ("OptionalProperty", python::no_init)
        .def("check", check<const Property&>)
        .def("get", get<const Property&>,
            python::return_value_policy<python::reference_existing_object>());

    python::class_<optional<const Residue&>> ("OptionalResidue", python::no_init)
        .def("check", check<const Residue&>)
        .def("get", get<const Residue&>,
            python::return_value_policy<python::reference_existing_object>());

    /**************************************************************************
     * Property
     **************************************************************************/
    python::class_<Property> ("Property", python::no_init)
        .def(python::init<bool>())
        .def(python::init<std::string>())
        .def(python::init<Vector3D>())
        .def(python::init<int>())
        .def(python::init<double>())
        .def("kind", &Property::kind)
        .def("as_bool", &Property::as_bool)
        .def("as_double", &Property::as_double)
        .def("as_vector3d", &Property::as_vector3d)
        .def("as_string", &Property::as_string,
            python::return_value_policy<python::copy_const_reference>());

    python::enum_<Property::Kind>("Kind")
        .value("BOOL", Property::BOOL)
        .value("DOUBLE", Property::DOUBLE)
        .value("STRING", Property::STRING)
        .value("VECTOR3D", Property::VECTOR3D);

    /**************************************************************************
     * Vector-likes
     **************************************************************************/
    python::class_<Bond, noncopyable>("Bond", python::no_init)
        .def("__getitem__", &get_index<Bond,size_t,2>);

    python::class_<std::vector<Bond>>("BondVec")
        .def(python::vector_indexing_suite<std::vector<Bond> >());

    python::class_<Angle, noncopyable>("Angle", python::no_init)
        .def("__getitem__", &get_index<Angle,size_t,3>);

    python::class_<std::vector<Angle>>("AngleVec")
        .def(python::vector_indexing_suite<std::vector<Angle> >());

    python::class_<Dihedral, noncopyable>("Dihedral", python::no_init)
        .def("__getitem__", &get_index<Dihedral,size_t,4>);

    python::class_<std::vector<Dihedral>>("DihedralVec")
        .def(python::vector_indexing_suite<std::vector<Dihedral> >());

    python::class_<Improper, noncopyable>("Improper", python::no_init)
        .def("__getitem__", &get_index<Improper,size_t,4>);

    python::class_<std::vector<Improper>>("ImproperVec")
        .def(python::vector_indexing_suite<std::vector<Improper> >());

    python::class_<Vector3D>("Vector3D")
        .def(python::init<double,double,double>())
        .def("__getitem__", &get_index<Vector3D,double,3>)
        .def("norm", &Vector3D::norm);

    python::class_<std::vector<Vector3D>>("PositionVec")
        .def(python::vector_indexing_suite<std::vector<Vector3D> >());

    /**************************************************************************
     * Residue
     **************************************************************************/
    python::class_<Residue, noncopyable>("Residue", python::init<std::string>())
        .def(python::init<std::string, int>())
        .def("size", &Residue::size)
        .def("name", &Residue::name,
            python::return_value_policy<python::copy_const_reference>())
        .def("atoms", python::range(&Residue::cbegin,
                                    &Residue::cend))
        .def("contains", &Residue::contains)
        .def("id", &Residue::id);
    python::def("get", getp<Residue>);

    /**************************************************************************
     * Topology
     **************************************************************************/
    python::class_<Topology, noncopyable>("Topology")
        .def("residue", &Topology::residue,
            python::return_value_policy<python::reference_existing_object>())
        .def("are_linked", &Topology::are_linked)
        .def("bonds", &Topology::bonds,
            python::return_value_policy<python::reference_existing_object>())
        .def("angles", &Topology::angles,
            python::return_value_policy<python::reference_existing_object>())
        .def("dihedrals", &Topology::dihedrals,
            python::return_value_policy<python::reference_existing_object>())
        .def("impropers", &Topology::impropers,
            python::return_value_policy<python::reference_existing_object>());

    /**************************************************************************
     * Frame
     **************************************************************************/
    python::class_<Frame, noncopyable>("Frame")
        .def("size", &Frame::size)
        .def("atoms", python::range(&Frame::cbegin,
                                    &Frame::cend))
        .def("topology", &Frame::topology,
            python::return_value_policy<python::reference_existing_object>())
        .def("distance", &Frame::distance)
        .def("angle", &Frame::angle)
        .def("dihedral", &Frame::dihedral)
        .def("out_of_plane", &Frame::out_of_plane);
    python::def("get", getp<Frame>);

    /**************************************************************************
     * Atom
     **************************************************************************/
    python::class_<Atom>("Atom", python::no_init)
        .def("name", &Atom::name,
            python::return_value_policy<python::copy_const_reference>())
        .def("type", &Atom::type,
            python::return_value_policy<python::copy_const_reference>())
        .def("mass", &Atom::mass)
        .def("charge", &Atom::charge)
        .def("full_name", &Atom::full_name)
        .def("vdw_radius", &Atom::vdw_radius)
        .def("covalent_radius", &Atom::covalent_radius)
        .def("atomic_number", &Atom::atomic_number);
    python::def("get", getp<Atom>);
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
    if (PyImport_AppendInittab("lemon", initlemon) == -1) {
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
    python::object result = python::exec(str.c_str(), global, global);

    // Obtained derived class from python
    python::object PythonDerived = global[py_derive.c_str()];
    python::object py_base = PythonDerived();
    LemonPythonBase& py = python::extract<LemonPythonBase&>(py_base);

    auto worker = [&py](chemfiles::Frame complex, const std::string& pdbid) {
        try {
            return py.worker(complex);
        } catch (...) {
            PyObject *ptype, *pvalue, *ptraceback;
            PyErr_Fetch(&ptype, &pvalue, &ptraceback);
            return pdbid + " " + PyString_AsString(pvalue) + "\n";
        }
    };

    lemon::launch<lemon::print_combine>(o, worker, std::cout);

    return 0;
}
