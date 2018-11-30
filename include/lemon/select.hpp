#ifndef LEMON_SELECT_HPP
#define LEMON_SELECT_HPP

#include <set>

#include "chemfiles/Frame.hpp"
#include "lemon/residue_name.hpp"

namespace lemon {

const ResidueNameSet common_peptides{
    {"CSD"}, {"PCA"}, {"DLE"}, {"KCX"}, {"CAS"}, {"CSO"}, {"PTR"},
    {"CME"}, {"SAH"}, {"TPO"}, {"SEP"}, {"MLY"}, {"HYP"}, {"MSE"},
    {"CYS"}, {"TRP"}, {"MET"}, {"HIS"}, {"TYR"}, {"GLN"}, {"PHE"},
    {"ASN"}, {"PRO"}, {"ARG"}, {"THR"}, {"ASP"}, {"ILE"}, {"LYS"},
    {"SER"}, {"GLU"}, {"VAL"}, {"GLY"}, {"ALA"}, {"LEU"},
};

/*!
 *  \addtogroup select
 *  @{
 */

//! \brief Select residues by removing them based on a criterion
namespace select {

inline std::set<size_t> small_molecules(const chemfiles::Frame& input,
                                        size_t min_atoms = 10) {

    const auto& residues = input.topology().residues();
    std::set<size_t> selected_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const auto& residue = residues[selected_residue];

        // Quick check, filters out many ions and small molecules before more
        // expensive checks.
        // For example, water
        if (residue.size() < min_atoms) {
            continue;
        }

        const auto& composition_type =
            residue.get("composition_type")->as_string();
        if (composition_type != "NON-POLYMER" && composition_type != "OTHER" &&
            composition_type != "PEPTIDE-LIKE") {
            continue;
        }

        size_t num_heavy_atoms = std::count_if(
            std::begin(residue), std::end(residue), [&input](size_t index) {
                return *(input[index].atomic_number()) != 1;
            });

        if (num_heavy_atoms < min_atoms) {
            continue;
        }

        selected_residues.insert(selected_residue);
    }

    return selected_residues;
}

inline std::set<size_t> metal_ions(const chemfiles::Frame& input) {

    const auto& residues = input.topology().residues();

    std::set<size_t> selected_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const auto& residue = residues[selected_residue];

        if (residue.size() == 1 && input[*residue.begin()].charge() > 0.0) {
            selected_residues.insert(selected_residue);
        }
    }

    return selected_residues;
}

inline std::set<size_t> nucleic_acids(const chemfiles::Frame& input) {

    const auto& residues = input.topology().residues();

    std::set<size_t> selected_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const chemfiles::Residue& residue = residues[selected_residue];
        const auto& comp_type = residue.get("composition_type")->as_string();

        if (comp_type.find("DNA") == std::string::npos &&
            comp_type.find("RNA") == std::string::npos) {
            continue;
        }

        selected_residues.insert(selected_residue);
    }

    return selected_residues;
}

inline std::set<size_t> peptides(const chemfiles::Frame& input) {
    const auto& residues = input.topology().residues();

    std::set<size_t> selected_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const chemfiles::Residue& residue = residues[selected_residue];
        const auto& comp_type = residue.get("composition_type")->as_string();

        if (comp_type.find("PEPTIDE") == std::string::npos ||
            comp_type == "PEPTIDE-LIKE") {
            continue;
        }

        selected_residues.insert(selected_residue);
    }

    return selected_residues;
}

inline std::set<size_t> specific_residues(
    const chemfiles::Frame& input, const ResidueNameSet& resnames) {

    const auto& residues = input.topology().residues();
    std::set<size_t> selected_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const auto& residue = residues[selected_residue];

        if (resnames.count({residue.name()}) == 0) {
            continue;
        }

        selected_residues.insert(selected_residue);
    }

    return selected_residues;
}

} // namespace select
/*! @} End of Doxygen Groups*/

} // namespace lemon

#endif
