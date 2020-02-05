#ifndef LEMON_SELECT_HPP
#define LEMON_SELECT_HPP

#include <list>
#include <set>
#include <unordered_set>

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#include <chemfiles/Frame.hpp>
LEMON_EXTERNAL_FILE_POP

#include "lemon/constants.hpp"
#include "lemon/residue_name.hpp"

namespace lemon {

//! Functions to select various residue based on a given criterion
namespace select {

//! Select small molecules in a given frame
//!
//! Use this function to find small molecules in a given `frame`. A small
//! molecule is defined as an entity that has a given [chemical composition]
//! (http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.type.html).
//! Also, the selected entity must have a specified number of atoms (default
//! 10), so that common residues such as water and metal ions are not selected.
//! \param [in] frame The entry containing molecules of interest.
//! \param [out] output A container filled with the small molecule residue ids
//! \param [in] types A set of `std::string` containing the accepted chemical
//!  chemical composition. Defaults are *NON-POLYMER*, *OTHER*, *PEPTIDE-LIKE*
//! \param [in] min_heavy_atoms The minimum number of non-hydrogen atoms for a
//!  residue to be classed as a small molecule.
//! \return The number of added residue IDs
template <typename Container>
inline size_t small_molecules(
    const chemfiles::Frame& frame, Container& output,
    const std::unordered_set<std::string>& types = small_molecule_types,
    size_t min_heavy_atoms = 10) {

    using chemfiles::Property;

    const auto& residues = frame.topology().residues();
    auto initialize_size = output.size();

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const auto& residue = residues[selected_residue];

        // Quick check, filters out many ions and small molecules before more
        // expensive checks.
        // For example, water
        if (residue.size() < min_heavy_atoms) {
            continue;
        }

        auto composition_type =
            residue.get<Property::STRING>("composition_type").value_or("");

        if (!types.count(composition_type)) {
            continue;
        }

        auto num_heavy_atoms = std::count_if(
            std::begin(residue), std::end(residue), [&frame](size_t index) {
                return *(frame[index].atomic_number()) != 1;
            });

        if (static_cast<size_t>(num_heavy_atoms) < min_heavy_atoms) {
            continue;
        }

        output.insert(output.end(), selected_residue);
    }

    return output.size() - initialize_size;
}

//! Select small molecules in a given frame
//!
//! Use this function to find small molecules in a given `frame`. A small
//! molecule is defined as an entity that has a given [chemical composition]
//! (http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.type.html).
//! Also, the selected entity must have a specified number of atoms (default
//! 10), so that common residues such as water and metal ions are not selected.
//! \param [in] frame The entry containing molecules of interest.
//! \param [in] types A set of `std::string` containing the accepted chemical
//!  chemical composition. Defaults are *NON-POLYMER*, *OTHER*, *PEPTIDE-LIKE*
//! \param [in] min_heavy_atoms The minimum number of non-hydrogen atoms for a
//!  residue to be classed as a small molecule.
//! \return A container filled with the small molecule residue ids
template <typename Container = std::list<size_t>>
inline Container small_molecules(
    const chemfiles::Frame& frame,
    const std::unordered_set<std::string>& types = small_molecule_types,
    size_t min_heavy_atoms = 10) {
    Container contain;
    small_molecules(frame, contain, types, min_heavy_atoms);
    return contain;
}

//! Select metal ions in a given frame
//!
//! This function populates the residue IDs of metal ions. We define a metal ion
//! as a residue with a single, positively charged ion.
//! \param [in] frame The entry containing metal ions of interest.
//! \param [out] output A container of residue IDs which are metals.
//! \return The number of added residue IDs
template <typename Container>
inline size_t metal_ions(const chemfiles::Frame& frame, Container& output) {

    const auto& residues = frame.topology().residues();
    auto initialize_size = output.size();

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const auto& residue = residues[selected_residue];

        if (residue.size() == 1 && frame[*residue.begin()].charge() > 0.0) {
            output.insert(output.end(), selected_residue);
        }
    }

    return output.size() - initialize_size;
}

//! Select metal ions in a given frame
//!
//! This function returns the residue IDs of metal ions. We define a metal ion
//! as a residue with a single, positively charged ion.
//! \param [in] frame The entry containing metal ions of interest.
//! \return output A container of residue IDs which are metals.
template <typename Container = std::list<size_t>>
inline Container metal_ions(const chemfiles::Frame& frame) {
    Container contain;
    metal_ions(frame, contain);
    return contain;
}

//! Select nucleic acid residues in a given frame
//!
//! This function populates the residue IDs of nucleic acid residues.
//! We define a nucleic acid as a residue with a chemical composition
//! containing the *RNA* or *DNA* substring.
//! \param [in] frame The entry containing nucleic acid residues.
//! \param [out] output A container of residue IDs which have nucleic acid
//! linkage \return The number of added residue IDs
template <typename Container>
inline size_t nucleic_acids(const chemfiles::Frame& frame, Container& output) {

    using chemfiles::Property;

    const auto& residues = frame.topology().residues();
    auto initialize_size = output.size();

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const chemfiles::Residue& residue = residues[selected_residue];
        auto comp_type =
            residue.get<Property::STRING>("composition_type").value_or("");

        if (comp_type.find("DNA") == std::string::npos &&
            comp_type.find("RNA") == std::string::npos) {
            continue;
        }

        output.insert(output.end(), selected_residue);
    }

    return output.size() - initialize_size;
}

//! Select nucleic acid residues in a given frame
//!
//! This function returns the residue IDs of nucleic acid residues.
//! We define a nucleic acid as a residue with a chemical composition
//! containing the *RNA* or *DNA* substring.
//! \param [in] frame The entry containing nucleic acid residues.
//! \return output A container of residue IDs which have nucleic acid linkage
template <typename Container = std::list<size_t>>
inline Container nucleic_acids(const chemfiles::Frame& frame) {
    Container contain;
    nucleic_acids(frame, contain);
    return contain;
}

//! Select peptide residues in a given frame
//!
//! This function populates the residue IDs of peptide residues.
//! We define a peptided as a residue with a chemical composition
//! containing the *PEPTIDE* substring which is not *PEPTIDE-LIKE*.
//! \param [in] frame The entry containing peptide residues.
//! \param [out] output A container of residue IDs which have peptide linkage
//! \return The number of added residue IDs
template <typename Container>
inline size_t peptides(const chemfiles::Frame& frame, Container& output) {

    using chemfiles::Property;

    const auto& residues = frame.topology().residues();
    auto initialize_size = output.size();

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const chemfiles::Residue& residue = residues[selected_residue];
        auto comp_type =
            residue.get<Property::STRING>("composition_type").value_or("");

        if (comp_type.find("PEPTIDE") == std::string::npos ||
            comp_type == "PEPTIDE-LIKE") {
            continue;
        }

        output.insert(output.end(), selected_residue);
    }

    return output.size() - initialize_size;
}

//! Select peptide residues in a given frame
//!
//! This function returns the residue IDs of peptide residues.
//! We define a peptided as a residue with a chemical composition
//! containing the *PEPTIDE* substring which is not *PEPTIDE-LIKE*.
//! \param [in] frame The entry containing peptide residues.
//! \return A container of residue IDs which have peptide linkage
template <typename Container = std::list<size_t>>
inline Container peptides(const chemfiles::Frame& frame) {
    Container contain;
    peptides(frame, contain);
    return contain;
}

//! Select residues with a given name in a given frame
//!
//! This function populates the residue IDs of peptides matching a given name
//! set. \param [in] frame The entry containing residues of interest. \param
//! [out] output The residue IDs with names that match those in `resnames`
//! \param [in] resnames The set of residue names of interest.
//! \return The number of added residue IDs
template <typename Container>
inline size_t specific_residues(const chemfiles::Frame& frame,
                                Container& output,
                                const ResidueNameSet& resnames) {

    const auto& residues = frame.topology().residues();
    auto initialize_size = output.size();

    for (size_t selected_residue = 0; selected_residue < residues.size();
         ++selected_residue) {
        const auto& residue = residues[selected_residue];

        if (resnames.count({residue.name()}) == 0) {
            continue;
        }

        output.insert(output.end(), selected_residue);
    }

    return output.size() - initialize_size;
}

//! Select residues with a given name in a given frame
//!
//! This function returns the residue IDs of peptides matching a given name set.
//! \param [in] frame The entry containing residues of interest.
//! \param [in] resnames The set of residue names of interest.
//! \return The residue IDs with names that match those in `resnames`
template <typename Container = std::list<size_t>>
inline Container specific_residues(const chemfiles::Frame& frame,
                                   const ResidueNameSet& resnames) {
    Container contain;
    specific_residues(frame, contain, resnames);
    return contain;
}

} // namespace select

} // namespace lemon

#endif
