#include <string>
#include <iostream>

#include "lemon/lemon.hpp"
#include "lemon/geometry.hpp"
#include "lemon/tmalign.hpp"
#include "lemon/xscore.hpp"

#include "chemfiles/File.hpp"

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

namespace lemon {

struct LemonPythonWrap : LemonPythonBase {
    virtual std::string worker(const chemfiles::Frame* frame,
                               const std::string& pdbid) override {
        py::gil_scoped_release release;
        {
            py::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(
                std::string,
                LemonPythonBase,
                worker,
                frame,
                pdbid
            );
        }
    }

    virtual void finalize() override {
        PYBIND11_OVERLOAD(
            void,
            LemonPythonBase,
            finalize,
        );
    }
};

template<typename T, typename ret, int maxval>
ret get_index(const T& v, int i) {
    if (i >= 0) {
        if (i >= maxval) {
            throw pybind11::index_error("index is too large");
        }
        return v[static_cast<size_t>(i)];
    }
    if (-i > maxval) {
        throw pybind11::index_error("index is too small");
    }
    return v[static_cast<size_t>(maxval + i)];
}

template<typename T, typename ret>
ret get_index(const T& v, int i) {
    if (i >= 0) {
        if (i >= static_cast<int>(v.size())) {
            throw pybind11::index_error("index is too large");
        }
        return v[static_cast<size_t>(i)];
    }
    if (-i > static_cast<int>(v.size())) {
        throw pybind11::index_error("index is too small");
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
        throw pybind11::index_error("Cannot dereference nullopt");
    }
    chemfiles::unreachable();
}

void translate(lemon::geometry::geometry_error const& e)
{
    // Use the Python 'C' API to set up an exception object
    auto msg = std::string("Geometry Error: ") + e.what();
    PyErr_SetString(PyExc_RuntimeError, msg.c_str());
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

template<typename T>
std::string to_string(const T& v) {
    std::stringstream ss;
    ss << v;
    return ss.str();
}

chemfiles::Frame open_model_in_file(const std::string& filename, size_t index) {
    chemfiles::Trajectory traj(filename);
    return traj.read_step(index);
}

chemfiles::Frame* open_file(const std::string& filename) {
    chemfiles::Trajectory traj(filename);
    return new chemfiles::Frame(traj.read());
}

void write_file(const chemfiles::Frame& frame, const std::string& filename) {
    chemfiles::Trajectory traj(filename, chemfiles::File::Mode::WRITE);
    traj.write(frame);
}

void append_file(const chemfiles::Frame& frame, const std::string& filename) {
    chemfiles::Trajectory traj(filename, chemfiles::File::Mode::APPEND);
    traj.write(frame);
}

void run_lemon_workflow(LemonPythonBase& py, const std::string& p, size_t threads) {
    auto worker = [&py](chemfiles::Frame complex, const std::string& pdbid) {
        try {
            return py.worker(&complex, pdbid);
        } catch (py::error_already_set& err) {
            return pdbid + " " + err.what() + "\n";
        } catch (py::cast_error& err) {
            return pdbid + " Problem with type: " + err.what() + "\n";
        } catch (std::exception& err) {
            return pdbid + " " + err.what() + "\n";
        } catch (...) {
            return pdbid + " unknown error." + "\n";
        }
    };

    print_combine combiner(std::cout);
    lemon::run_parallel(worker, p, combiner, threads);
}
}

// Pack the Base class wrapper into a module
PYBIND11_MODULE(lemon, m) {
    using namespace chemfiles;
    using namespace lemon;

    py::class_<LemonPythonBase, LemonPythonWrap>(m, "Workflow")
        .def(py::init<>())
        .def("worker", &LemonPythonBase::worker)
        .def("finalize", &LemonPythonBase::finalize);

    m.def("launch", run_lemon_workflow);

    /**************************************************************************
     * Optional
     **************************************************************************/
    py::class_<optional<uint64_t>> (m, "OptionalUInt64")
        .def("check", check<uint64_t>)
        .def("get", get<uint64_t>);

    py::class_<optional<std::string>> (m, "OptionalString")
        .def("check", check<std::string>)
        .def("get", get<std::string>);

    py::class_<optional<double>> (m, "OptionalDouble")
        .def("check", check<double>)
        .def("get", get<double>);

    py::class_<optional<const Property&>> (m, "OptionalProperty")
        .def("check", check<const Property&>)
        .def("get", get<const Property&>,
            py::return_value_policy::reference_internal);

    py::class_<optional<const Residue&>> (m, "OptionalResidue")
        .def("check", check<const Residue&>)
        .def("get", get<const Residue&>,
            py::return_value_policy::reference_internal);

    /**************************************************************************
     * Property
     **************************************************************************/
    py::class_<Property> (m, "Property")
        .def(py::init<bool>())
        .def(py::init<std::string>())
        .def(py::init<Vector3D>())
        .def(py::init<int>())
        .def(py::init<double>())
        .def("kind", &Property::kind)
        .def("as_bool", &Property::as_bool)
        .def("as_double", &Property::as_double)
        .def("as_vector3d", &Property::as_vector3d)
        .def("as_string", &Property::as_string);

    py::enum_<Property::Kind>(m, "Kind")
        .value("BOOL", Property::BOOL)
        .value("DOUBLE", Property::DOUBLE)
        .value("STRING", Property::STRING)
        .value("VECTOR3D", Property::VECTOR3D);

    /**************************************************************************
     * Vector-likes
     **************************************************************************/
    py::class_<Bond>(m, "Bond")
        .def("__getitem__", &get_index<Bond,size_t,2>);

    py::bind_vector<std::vector<Bond>>(m, "BondVec");

    py::class_<Angle>(m, "Angle")
        .def("__getitem__", &get_index<Angle,size_t,3>);

    py::bind_vector<std::vector<Angle>>(m, "AngleVec");

    py::class_<Dihedral>(m, "Dihedral")
        .def("__getitem__", &get_index<Dihedral,size_t,4>);

    py::bind_vector<std::vector<Dihedral>>(m, "DihedralVec");

    py::class_<Improper>(m, "Improper")
        .def("__getitem__", &get_index<Improper,size_t,4>);

    py::bind_vector<std::vector<Improper>>(m, "ImproperVec");

    py::class_<Vector3D>(m, "Vector3D")
        .def(py::init<double,double,double>())
        .def("__getitem__", &get_index<Vector3D,double,3>)
        .def("norm", &Vector3D::norm);

    py::bind_vector<std::vector<Vector3D>>(m, "PositionVec");

    /**************************************************************************
     * Residue
     **************************************************************************/
    chemfiles::optional<const chemfiles::Property&> (Residue::*residue_get)
        (const std::string&) const = &Residue::get;

    py::class_<Residue>(m, "Residue")
        .def(py::init<std::string>())
        .def(py::init<std::string, size_t>())
        .def("__len__", &Residue::size)
        .def("__iter__", [](const Residue& v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("name", &Residue::name)
        .def("contains", &Residue::contains)
        .def("id", &Residue::id)
        .def("get", residue_get);

    py::class_<std::vector<Residue>>(m, "ResidueVec");

    /**************************************************************************
     * Topology
     **************************************************************************/
    py::class_<Topology>(m, "Topology")
        .def("__len__", &Topology::size)
        .def("__get_item__", get_index<Topology,const Atom&>,
            py::return_value_policy::reference_internal)
        .def("__iter__", [](const Topology& v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("residue", &Topology::residue,
            py::return_value_policy::reference_internal)
        .def("residues", &Topology::residues,
            py::return_value_policy::reference_internal)
        .def("residue_for_atom", &Topology::residue_for_atom)
        .def("are_linked", &Topology::are_linked)
        .def("bonds", &Topology::bonds,
            py::return_value_policy::reference_internal)
        .def("angles", &Topology::angles,
            py::return_value_policy::reference_internal)
        .def("dihedrals", &Topology::dihedrals,
            py::return_value_policy::reference_internal)
        .def("impropers", &Topology::impropers,
            py::return_value_policy::reference_internal);

    /**************************************************************************
     * Frame
     **************************************************************************/
    chemfiles::optional<const chemfiles::Property&> (Frame::*frame_get)
        (const std::string&) const = &Frame::get;

    py::class_<Frame>(m, "Frame")
        .def("__len__", &Frame::size)
        .def("__get_item__", get_index<Frame,const Atom&>,
            py::return_value_policy::reference_internal)
        .def("__iter__", [](const Frame& v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("topology", &Frame::topology,
            py::return_value_policy::reference_internal)
        .def("distance", &Frame::distance)
        .def("angle", &Frame::angle)
        .def("dihedral", &Frame::dihedral)
        .def("out_of_plane", &Frame::out_of_plane)
        .def("get", frame_get);

    /**************************************************************************
     * Atom
     **************************************************************************/
    chemfiles::optional<const chemfiles::Property&> (Atom::*atom_get)
        (const std::string&) const = &Atom::get;

    py::class_<Atom>(m, "Atom")
        .def("name", &Atom::name)
        .def("type", &Atom::type)
        .def("mass", &Atom::mass)
        .def("charge", &Atom::charge)
        .def("full_name", &Atom::full_name)
        .def("vdw_radius", &Atom::vdw_radius)
        .def("covalent_radius", &Atom::covalent_radius)
        .def("atomic_number", &Atom::atomic_number)
        .def("get", atom_get);

    /**************************************************************************
     * Residue Name
     **************************************************************************/
    py::class_<ResidueName>(m, "ResidueName")
        .def(py::init<const std::string&>())
        .def("__str__",[](const ResidueName& v){
            return to_string(v);
        })
        .def("__repl__",[](const ResidueName& v){
            return "<ResidueName {" + to_string(v) + "}";
        });

    typedef std::pair<ResidueNameSet::iterator, bool> rns_insert_ret;
    py::class_<rns_insert_ret>(m, "ResidueNameRet");

    rns_insert_ret (ResidueNameSet::*rns_insert)
        (const ResidueNameSet::value_type&) = &ResidueNameSet::insert;
    py::class_<ResidueNameSet>(m, "ResidueNameSet")
        .def("__len__", &ResidueNameSet::size)
        .def("__iter__", [](const ResidueNameSet& v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("append", rns_insert)
        .def("__str__",[](const ResidueNameSet& v){
            return to_string(v);
        })
        .def("__repl__",[](const ResidueNameSet& v){
            return "ResidueNameSet {" + to_string(v) + "}";
        });

    py::bind_map<ResidueNameCount>(m, "ResidueNameCount");

    void (default_id_list::*push_back)(const default_id_list::value_type&) =
        &default_id_list::push_back;

    py::class_<default_id_list>(m, "ResidueIDs")
        .def("__iter__", [](const default_id_list& v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("append", push_back)
        .def("__str__",[](const default_id_list& v){
            return to_string(v);
        })
        .def("__repl__",[](const default_id_list& v){
            return "ResidueIDs {" + to_string(v) + "}";
        });

    /**************************************************************************
     * Constants
     **************************************************************************/
    py::class_<std::unordered_set<std::string>>(m, "StringSet");
    m.attr("small_molecule_types") = small_molecule_types;

    m.attr("common_peptides") = common_peptides;
    m.attr("common_cofactors") = common_cofactors;
    m.attr("common_fatty_acids") = common_fatty_acids;
    m.attr("proline_res") = proline_res;

    /**************************************************************************
     * Select
     **************************************************************************/

    // Returns a new object    
    default_id_list (*small_molecules)(const Frame&,
                                       const std::unordered_set<std::string>&,
                                       size_t) =
        &select::small_molecules;

    m.def("select_small_molecules", small_molecules);

    default_id_list (*metal_ions)(const Frame&) = &select::metal_ions;
    m.def("select_metal_ions", metal_ions);

    default_id_list (*nucleic_acids)(const Frame&) = &select::nucleic_acids;
    m.def("select_nucleic_acids", nucleic_acids);

    default_id_list (*peptides)(const Frame&) = &select::peptides;
    m.def("select_peptides", peptides);

    default_id_list (*specific_residues)(const Frame&, const ResidueNameSet&) =
        &select::specific_residues;
    m.def("select_specific_residues", specific_residues);

    // Inplace
    size_t (*small_molecules_i)(const Frame&, default_id_list&,
                                const std::unordered_set<std::string>&,
                                size_t) =
        &select::small_molecules;

    m.def("select_small_molecules", small_molecules_i);

    size_t (*metal_ions_i)(const Frame&, default_id_list&) =
         &select::metal_ions;
    m.def("select_metal_ions", metal_ions_i);

    size_t (*nucleic_acids_i)(const Frame&, default_id_list&) =
        &select::nucleic_acids;
    m.def("select_nucleic_acids", nucleic_acids_i);

    size_t (*peptides_i)(const Frame&, default_id_list&) =
        &select::peptides;
    m.def("select_peptides", peptides_i);

    size_t (*specific_residues_i)(const Frame&, default_id_list&,
                                  const ResidueNameSet&) =
        &select::specific_residues;
    m.def("select_specific_residues", specific_residues_i);

    /**************************************************************************
     * Count
     **************************************************************************/
    m.def("count_atomic_property", count::atom_property);
    m.def("count_residue_property", count::residue_property);
    m.def("count_print_residue_names",
        count::print_residue_names<default_id_list>);

    void (*residues1)(const Frame&, ResidueNameCount&) = &count::residues;
    m.def("count_residues", residues1);

    void (*residues2)(const Frame&, const default_id_list&, ResidueNameCount&) =
        &count::residues;
    m.def("count_residues", residues2);

    /**************************************************************************
     * Prune
     **************************************************************************/
    m.def("prune_identical_residues",
        prune::identical_residues<default_id_list>);
    m.def("prune_cofactors", prune::cofactors<default_id_list>);
    m.def("keep_interactions", prune::keep_interactions<default_id_list>);
    m.def("remove_interactions", prune::remove_interactions<default_id_list>);

    /**************************************************************************
     * Separate
     **************************************************************************/
    m.def("separate_residues", separate::residues<default_id_list>);
    m.def("separate_protein_and_ligand", separate::protein_and_ligand);

    /**************************************************************************
     * geometry
     **************************************************************************/
    m.def("protein_bond_name", geometry::protein::bond_name);
    m.def("protein_angle_name", geometry::protein::angle_name);
    m.def("protein_dihedral_name", geometry::protein::dihedral_name);
    m.def("protein_improper_name", geometry::protein::improper_name);
    py::register_exception<geometry::geometry_error>(m,"GeometryError");

    /**************************************************************************
     * Vina Score
     **************************************************************************/
    py::class_<xscore::VinaScore>(m,"VinaScore")
        .def_readonly("g1", &xscore::VinaScore::g1)
        .def_readonly("g2", &xscore::VinaScore::g2)
        .def_readonly("rep", &xscore::VinaScore::rep)
        .def_readonly("hydrophobic", &xscore::VinaScore::hydrophobic)
        .def_readonly("hydrogen", &xscore::VinaScore::hydrogen);

    m.def("vina_score", xscore::vina_score<default_id_list>);

    /**************************************************************************
     * TMAlign
     **************************************************************************/
    py::class_<tmalign::TMResult>(m,"TMResult")
        .def_readonly("score", &tmalign::TMResult::score)
        .def_readonly("rmsd", &tmalign::TMResult::rmsd)
        .def_readonly("aligned", &tmalign::TMResult::aligned);

    m.def("TMscore", tmalign::TMscore);

    /**************************************************************************
    * File IO
    ***************************************************************************/
    m.def("open_model_in_file", open_model_in_file);
    m.def("open_file", open_file);
    m.def("write_file", write_file);
    m.def("append_file", append_file);
}

#pragma clang diagnostic pop
