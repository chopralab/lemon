#ifndef LEMON_SELECT_HPP
#define LEMON_SELECT_HPP

#include <set>
#include <unordered_set>

#include "chemfiles/Frame.hpp"
#include "lemon/residue_name.hpp"
#include "lemon/constants.hpp"

namespace lemon {

//! Functions to select various residue based on a given criterion
namespace select {

//! Select small molecules in a given frame
//!
//! Use this function to find small molecules in a given `frame`. A small
//! molecule is defined as an entity that has a given [chemical composition]
//! (http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.type.html).
//! Also, the selected entity must have a specified number of atoms (default 10),
//! so that common residues such as water and metal ions are not selected.
//! \param [in] frame The complex containing molecules of interest.
//! \param [in] types A set of `std::string` containing the accepted chemical
//!  chemical composition. Defaults are *NON-POLYMER*, *OTHER*, *PEPTIDE-LIKE*
//! \param [in] min_heavy_atoms The minimum number of non-hydrogen atoms for a
//!  residue to be classed as a small molecule.
//! \return A set of residue IDs which match the supplied criterion.
inline std::set<size_t> small_molecules(const chemfiles::Frame& frame,
                                        const std::unordered_set<std::string>& types = small_molecule_types,
                                        size_t min_heavy_atoms = 10) {

    const auto& residues = frame.topology().residues();
    std::set<size_t> selected_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const auto& residue = residues[selected_residue];

        // Quick check, filters out many ions and small molecules before more
        // expensive checks.
        // For example, water
        if (residue.size() < min_heavy_atoms) {
            continue;
        }

        const auto& composition_type =
            residue.get("composition_type")->as_string();
        if (!types.count(composition_type)) {
            continue;
        }

        size_t num_heavy_atoms = std::count_if(
            std::begin(residue), std::end(residue), [&frame](size_t index) {
                return *(frame[index].atomic_number()) != 1;
            });

        if (num_heavy_atoms < min_heavy_atoms) {
            continue;
        }

        selected_residues.insert(selected_residue);
    }

    return selected_residues;
}

//! Select metal ions in a given frame
//!
//! This function returns the residue IDs of metal ions. We define a metal ion
//! as a residue with a single, positively charged ion.
//! \param [in] frame The complex containing metal ions of interest.
//! \return A set of residue IDs which match the supplied criterion.
inline std::set<size_t> metal_ions(const chemfiles::Frame& frame) {

    const auto& residues = frame.topology().residues();

    std::set<size_t> selected_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const auto& residue = residues[selected_residue];

        if (residue.size() == 1 && frame[*residue.begin()].charge() > 0.0) {
            selected_residues.insert(selected_residue);
        }
    }

    return selected_residues;
}

//! Select nucleic acid residues in a given frame
//!
//! This function returns the residue IDs of nucleic acid residues.
//! We define a nuecleic acid as a residue with a chemical composition
//! containing the *RNA* or *DNA* substring.
//! \param [in] frame The complex containing nucleic acid residues.
//! \return A set of residue IDs which match the supplied criterion.
inline std::set<size_t> nucleic_acids(const chemfiles::Frame& frame) {

    const auto& residues = frame.topology().residues();

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

//! Select peptide residues in a given frame
//!
//! This function returns the residue IDs of peptide residues.
//! We define a peptided as a residue with a chemical composition
//! containing the *PEPTIDE* substring which is not *PEPTIDE-LIKE*.
//! \param [in] frame The complex containing peptide residues.
//! \return A set of residue IDs which match the supplied criterion.
inline std::set<size_t> peptides(const chemfiles::Frame& frame) {
    const auto& residues = frame.topology().residues();

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

//! Select residues with a given name in a given frame
//!
//! This function returns the residue IDs of peptides matching a given name set.
//! \param [in] frame The complex containing residues of interest.
//! \param [in] resnames The set of residue names of interest.
//! \return A set of residue IDs which match the supplied criterion.
inline std::set<size_t> specific_residues(
    const chemfiles::Frame& frame, const ResidueNameSet& resnames) {

    const auto& residues = frame.topology().residues();
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

} // namespace lemon

#endif
