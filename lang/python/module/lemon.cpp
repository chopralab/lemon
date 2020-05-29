#include <mutex>
#include <string>
#include <iostream>

#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"
#include "lemon/geometry.hpp"
#include "lemon/tmalign.hpp"
#include "lemon/xscore.hpp"
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
#include <chemfiles.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
LEMON_EXTERNAL_FILE_POP

namespace py = pybind11;

namespace lemon {

struct LemonPythonWrap : LemonPythonBase {
    virtual std::string worker(const chemfiles::Frame* frame,
                               const std::string& pdbid) override {
        py::gil_scoped_acquire acq;
        try {
            PYBIND11_OVERLOAD_PURE(
                std::string,
                LemonPythonBase,
                worker,
                frame,
                pdbid
            );
        } catch (py::error_already_set& err) {
            std::cerr << pdbid + " " + err.what() + "\n";
        } catch (py::cast_error& err) {
            std::cerr << pdbid + " Problem with type: " + err.what() + "\n";
        } catch (std::exception& err) {
            std::cerr << pdbid + " " + err.what() + "\n";
        } catch (...) {
            std::cerr << pdbid + " unknown error." + "\n";
        }

        return std::string("");
    }

    virtual void finalize() override {
        PYBIND11_OVERLOAD(
            void,
            LemonPythonBase,
            finalize,
        );
    }
};

void run_lemon_workflow(LemonPythonBase& py, const std::string& p, size_t threads, const Entries& entries) {
    py::gil_scoped_release release;

    auto worker = [&py](chemfiles::Frame entry, const std::string& pdbid) {
        return py.worker(&entry, pdbid);
    };

    print_combine combiner(std::cout);
    lemon::run_parallel(worker, p, combiner, threads, entries);
    py.finalize();
}

using default_id_list = std::list<uint64_t>;
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

std::string to_string(const Affine& affine) {
    std::stringstream ss;
    
    ss << "[ T( " << affine.T[0] << ", " << affine.T[1] << ", " << affine.T[2]
       << ")";

    ss << " R( [" << affine.R[0][0] << ", " << affine.R[0][1] << ", "
       << affine.R[0][2] << "]" ;

    ss << "; [" << affine.R[1][0] << ", " << affine.R[1][1] << ", "
       << affine.R[1][2] << "]";

    ss << "; [" << affine.R[2][0] << ", " << affine.R[2][1] << ", "
       << affine.R[2][2] << "] )";

    return ss.str();
}

void translate(lemon::geometry::geometry_error const& e) {
    // Use the Python 'C' API to set up an exception object
    auto msg = std::string("Geometry Error: ") + e.what();
    PyErr_SetString(PyExc_RuntimeError, msg.c_str());
}

}

void add_lemon_features(py::module& m) {

    using namespace chemfiles;
    using namespace lemon;

    py::class_<LemonPythonBase, LemonPythonWrap>(m, "Workflow")
        .def(py::init<>())
        .def("worker", &LemonPythonBase::worker)
        .def("finalize", &LemonPythonBase::finalize);

    m.def("launch", run_lemon_workflow);
    m.def("launch", [](LemonPythonBase& py, const std::string& p, size_t threads){
        run_lemon_workflow(py, p, threads, Entries());
    });

    /**************************************************************************
     * Residue Name
     **************************************************************************/
    py::class_<ResidueName>(m, "ResidueName")
        .def(py::init<const std::string&>())
        .def("__str__",[](const ResidueName& v){
            return to_string(v);
        })
        .def("__repl__",[](const ResidueName& v){
            return "<ResidueName {" + to_string(v) + "}>";
        });

    /**************************************************************************
     * Constants
     **************************************************************************/
    m.attr("small_molecule_types") = small_molecule_types;

    m.attr("common_peptides") = common_peptides;
    m.attr("common_cofactors") = common_cofactors;
    m.attr("common_fatty_acids") = common_fatty_acids;
    m.attr("proline_res") = proline_res;

    /**************************************************************************
     * Select
     **************************************************************************/

    m.def("select_small_molecules", &select::small_molecules<default_id_list>);
    m.def("select_metal_ions", &select::metal_ions<default_id_list>);
    m.def("select_nucleic_acids", &select::nucleic_acids<default_id_list>);
    m.def("select_peptides", &select::peptides<default_id_list>);
    m.def("select_residue_ids", &select::residue_ids<default_id_list>);
    m.def("select_specific_residues",
          &select::specific_residues<default_id_list>);
    m.def("select_residue_property",
          &select::residue_property<default_id_list>);

    /**************************************************************************
     * Count
     **************************************************************************/
    m.def("count_atomic_property", count::atom_property);
    m.def("count_residue_property", count::residue_property);
    m.def("count_print_residue_names",
        count::print_residue_names<default_id_list>);

    ResidueNameCount& (*residues1)(const Frame&, ResidueNameCount&) =
        &count::residues;
    m.def("count_residues", residues1);

    ResidueNameCount& (*residues2)(const Frame&, const default_id_list&,
                                   ResidueNameCount&) =
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
    m.def("intersection", prune::intersection<default_id_list>);
    m.def("has_property", prune::has_property<default_id_list>);

    /**************************************************************************
     * Separate
     **************************************************************************/
    m.def("separate_residues", separate::residues<default_id_list>);

    m.def("separate_residues", [](const chemfiles::Frame& input,
                                  const default_id_list& accepted_residues,
                                  chemfiles::Frame& new_frame) {
        separate::residues(input, accepted_residues, new_frame);
    });

    m.def("separate_protein_and_ligand", separate::protein_and_ligand);

    m.def("separate_protein_and_ligand", [](const chemfiles::Frame& input,
                                            size_t ligand_id,
                                            double pocket_size,
                                            chemfiles::Frame& protein,
                                            chemfiles::Frame& ligand) {
        separate::protein_and_ligand(input, ligand_id, pocket_size, protein, ligand);
    });

    m.def("separate_protein_and_ligands", separate::protein_and_ligands<default_id_list>);

    m.def("separate_protein_and_ligands", [](const chemfiles::Frame& input,
                                             const default_id_list& ligand_ids,
                                             double pocket_size,
                                             chemfiles::Frame& protein,
                                             chemfiles::Frame& ligand) {
        separate::protein_and_ligands(input, ligand_ids, pocket_size, protein, ligand);
    });

    /**************************************************************************
     * geometry
     **************************************************************************/
    m.def("protein_bond_name", geometry::protein::bond_name);
    m.def("protein_angle_name", geometry::protein::angle_name);
    m.def("protein_dihedral_name", geometry::protein::dihedral_name);
    m.def("protein_improper_name", geometry::protein::improper_name);
    py::register_exception<geometry::geometry_error>(m,"GeometryError");

    /**************************************************************************
     * Matrix
     **************************************************************************/
    py::class_<Affine>(m, "Affine") // Python cannot modify or read this
        .def("__str__", [](const Affine& v) { return to_string(v); })
        .def("__repl__", [](const Affine& v) {
            return "<Affine {" + to_string(v) + "}>";
        });

    Affine (*kabsch_f)(Coordinates&, Coordinates&, double) = &kabsch;
    m.def("kabsch", kabsch_f);

    void (*align_f)(span<Vector3D>&, const Affine&) = &align;
    m.def("align", align_f);

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
        .def_readonly("aligned", &tmalign::TMResult::aligned)
        .def_readonly("affine", &tmalign::TMResult::affine);

    m.def("TMscore", tmalign::TMscore);
}
