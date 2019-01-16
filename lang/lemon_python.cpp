#include <string>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "lemon/lemon.hpp"

namespace python = boost::python;

namespace lemon {

class LemonPythonBase : public boost::noncopyable {
public:
    virtual ~LemonPythonBase() {}
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

template<typename T, typename ret>
ret get_index(const T& v, int i) {
    if (i >= 0) {
        if (i >= static_cast<int>(v.size())) {
            PyErr_SetString(PyExc_IndexError, "index is too large");
            python::throw_error_already_set();
        }
        return v[static_cast<size_t>(i)];
    }
    if (-i > static_cast<int>(v.size())) {
        PyErr_SetString(PyExc_IndexError, "index is too small");
        python::throw_error_already_set();
    }
    return v[static_cast<size_t>(static_cast<int>(v.size()) + i)];
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

using default_id_list = std::list<size_t>;
inline std::ostream& operator<<(std::ostream& os, const default_id_list& idlist) {
    os << '[';
    for (auto i : idlist) {
        os << i << ' ';
    }
    os << ']';
    return os;
}
}

// Pack the Base class wrapper into a module
BOOST_PYTHON_MODULE(lemon) {
    using namespace chemfiles;
    using namespace lemon;
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
            python::return_internal_reference<>());

    python::class_<optional<const Residue&>> ("OptionalResidue", python::no_init)
        .def("check", check<const Residue&>)
        .def("get", get<const Residue&>,
            python::return_internal_reference<>());

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
            python::return_internal_reference<>())
        .def("residue_for_atom", &Topology::residue_for_atom)
        .def("are_linked", &Topology::are_linked)
        .def("bonds", &Topology::bonds,
            python::return_internal_reference<>())
        .def("angles", &Topology::angles,
            python::return_internal_reference<>())
        .def("dihedrals", &Topology::dihedrals,
            python::return_internal_reference<>())
        .def("impropers", &Topology::impropers,
            python::return_internal_reference<>());

    /**************************************************************************
     * Frame
     **************************************************************************/
    python::class_<Frame, noncopyable>("Frame")
        .def("size", &Frame::size)
        .def("atoms", python::range(&Frame::cbegin,
                                    &Frame::cend))
        .def("__get_item__", get_index<Frame,const Atom&>,
            python::return_internal_reference<>())
        .def("topology", &Frame::topology,
            python::return_internal_reference<>())
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

    /**************************************************************************
     * Residue Name
     **************************************************************************/
    python::class_<ResidueName>("ResidueName", python::init<const std::string&>())
        .def(python::self_ns::str(python::self));

    python::class_<ResidueNameSet>("ResidueNameSet")
        .def(python::self_ns::str(python::self));

    python::class_<ResidueNameCount>("ResidueNameCount")
        .def(python::self_ns::str(python::self))
        .def(python::self += python::self);

    python::class_<default_id_list>("ResidueIDs")
        //.def(python::self_ns::str(python::self))
        .def("__iter__", python::range(&default_id_list::cbegin,
                                       &default_id_list::cend))
        .def("size", &default_id_list::size);

    /**************************************************************************
     * Constants
     **************************************************************************/
    python::class_<std::unordered_set<std::string>>("StringSet");
    python::scope().attr("small_molecule_types") = small_molecule_types;

    python::scope().attr("common_peptides") = common_peptides;
    python::scope().attr("common_cofactors") = common_cofactors;
    python::scope().attr("common_fatty_acids") = common_fatty_acids;

    /**************************************************************************
     * Select
     **************************************************************************/

    // Returns a new object    
    default_id_list (*small_molecules)(const Frame&,
                                       const std::unordered_set<std::string>&,
                                       size_t) =
        &select::small_molecules;

    python::def("select_molecules", small_molecules);

    default_id_list (*metal_ions)(const Frame&) = &select::metal_ions;
    python::def("select_metal_ions", metal_ions);

    default_id_list (*nucleic_acids)(const Frame&) = &select::nucleic_acids;
    python::def("select_nucleic_acids", nucleic_acids);

    default_id_list (*peptides)(const Frame&) = &select::peptides;
    python::def("select_peptides", peptides);

    default_id_list (*specific_residues)(const Frame&, const ResidueNameSet&) =
        &select::specific_residues;
    python::def("select_specific_residues", specific_residues);

    // Inplace
    size_t (*small_molecules_i)(const Frame&, default_id_list&,
                                const std::unordered_set<std::string>&,
                                size_t) =
        &select::small_molecules;

    python::def("select_molecules", small_molecules_i);

    size_t (*metal_ions_i)(const Frame&, default_id_list&) =
         &select::metal_ions;
    python::def("select_metal_ions", metal_ions_i);

    size_t (*nucleic_acids_i)(const Frame&, default_id_list&) =
        &select::nucleic_acids;
    python::def("select_nucleic_acids", nucleic_acids_i);

    size_t (*peptides_i)(const Frame&, default_id_list&) =
        &select::peptides;
    python::def("select_peptides", peptides_i);

    size_t (*specific_residues_i)(const Frame&, default_id_list&,
                                  const ResidueNameSet&) =
        &select::specific_residues;
    python::def("select_specific_residues", specific_residues_i);

    /**************************************************************************
     * Count
     **************************************************************************/
    python::def("count_altloc", count::altloc);
    python::def("count_bioassemblies", count::bioassemblies);
    python::def("print_residue_name_counts",
        count::print_residue_name_counts<default_id_list>);

    void (*residues1)(const Frame&, ResidueNameCount&) = &count::residues;
    python::def("count_residues", residues1);

    void (*residues2)(const Frame&, const default_id_list&, ResidueNameCount&) =
        &count::residues;
    python::def("count_residues", residues2);

    /**************************************************************************
     * Prune
     **************************************************************************/
    python::def("prune_identical_residues",
        prune::identical_residues<default_id_list>);
    python::def("prune_cofactors", prune::cofactors<default_id_list>);
    python::def("keep_interactions", prune::keep_interactions<default_id_list>);
    python::def("remove_interactions", prune::remove_interactions<default_id_list>);

    /**************************************************************************
     * Separate
     **************************************************************************/
    python::def("separate_residues", separate::residues<default_id_list>);
    python::def("separate_protein_and_ligand", separate::protein_and_ligand);

    /**************************************************************************
     * Vina Score
     **************************************************************************/
    python::class_<xscore::VinaScore>("VinaScore", python::no_init)
        .def_readonly("g1", &xscore::VinaScore::g1)
        .def_readonly("g2", &xscore::VinaScore::g2)
        .def_readonly("rep", &xscore::VinaScore::rep)
        .def_readonly("hydrophobic", &xscore::VinaScore::hydrophobic)
        .def_readonly("hydrogen", &xscore::VinaScore::hydrogen);

    python::def("vina_score", xscore::vina_score<default_id_list>);

    /**************************************************************************
     * TMAlign
     **************************************************************************/
    python::class_<tmalign::TMResult>("TMResult", python::no_init)
        .def_readonly("score", &tmalign::TMResult::score)
        .def_readonly("rmsd", &tmalign::TMResult::rmsd)
        .def_readonly("aligned", &tmalign::TMResult::aligned);

    python::def("TMscore", tmalign::TMscore);
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
    #if PY_MAJOR_VERSION >= 3
    if (PyImport_AppendInittab("lemon", PyInit_lemon) == -1) {
        std::cerr << "Failed to embed lemon in to builtin modules" << std::endl;
        return 1;
    }
    #else
    if (PyImport_AppendInittab("lemon", initlemon) == -1) {
        std::cerr << "Failed to embed lemon in to builtin modules" << std::endl;
        return 1;
    }
    #endif

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
    lemon::LemonPythonBase& py = python::extract<lemon::LemonPythonBase&>(py_base);

    auto worker = [&py](chemfiles::Frame complex, const std::string& pdbid) {
        try {
            return py.worker(complex);
        } catch (...) {
            PyObject *ptype, *pvalue, *ptraceback;
            PyErr_Fetch(&ptype, &pvalue, &ptraceback);
            #if PY_MAJOR_VERSION >= 3
            PyObject* pyStr = PyUnicode_AsEncodedString(pvalue, "utf-8", "Error ~");
            return pdbid + " " + PyBytes_AS_STRING(pyStr) + "\n";
            #else
            return pdbid + " " + PyString_AsString(pvalue) + "\n";
            #endif
        }
    };

    lemon::launch<lemon::print_combine>(o, worker, std::cout);

    return 0;
}
