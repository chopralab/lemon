#ifndef LEMON_PRUNE_HPP
#define LEMON_PRUNE_HPP

#include <chemfiles.hpp>
#include "lemon/residue_name.hpp"

namespace lemon {

//! Prune selected residues by removing them based on a criterion
namespace prune {

//! Remove residues which are biologic copies of one another in a crystal 
//!
//! Many crystal structures in the PDB contain two identical copies of a
//! biological macromolecule. Since these copies are functionally identical,
//! some users wishing to only analyze a unique set of protein chains may want
//! to remove the identical residue copy. This function performs this operation
//! on a given set of residue ids by comparing all the residues in the `frame`'s
//! biological assemblies. If a residue in one assembly has the same ID as a
//! residue in a different assembly, then the copied residue is removed.
//! \param [in] frame The `frame` containing residues of interest.
//! \param [in,out] residue_ids The residue IDs to be pruned
inline void identical_residues(const chemfiles::Frame& frame,
                               std::set<size_t>& residue_ids) {
    auto& residues = frame.topology().residues();

    auto it = residue_ids.begin();
    while (it != residue_ids.end()) {
        auto current_id = *it;
        const auto& res_current = residues[current_id];

        auto checking = it;
        ++checking;
        while (checking != residue_ids.end()) {
            auto check_current = checking++;
            auto check_id = *check_current;
            const auto& res_check = residues[check_id];

            // This is faster than a string comparison
            if (res_current.size() != res_check.size()) {
                continue;
            }

            // Are these named the same?
            if (res_current.name() == res_check.name()) {
                // We use bioassemblies because they are standardized by the
                // protein databank
                auto bio_current =
                    res_current.get("assembly")
                        ? res_current.get("assembly")->as_string()
                        : "";
                auto bio_check = res_check.get("assembly")
                                     ? res_check.get("assembly")->as_string()
                                     : "";

                if (bio_current != bio_check) {
                    residue_ids.erase(check_current);
                }
            }
        }

        ++it;
    }
}

//! Remove residues which are typically present in many crystal structures
//!
//! There are a common set of cofactors present in many crystal structures such
//! as sugars and fatty acids used to induce crystalization. As a result, some
//! users may remove these cofactors as they may match other criteria (such as
//! being a small molecule) set by the user. 
//! \param [in] frame The `frame` containing residues of interest.
//! \param [in,out] residue_ids The residue IDs to be pruned.
//! \param [in] rns The residue names that one wishes to remove from residue_ids.
inline void cofactors(const chemfiles::Frame& frame,
                      std::set<size_t>& residue_ids,
                      const ResidueNameSet& rns) {
    auto& residues = frame.topology().residues();

    auto it = residue_ids.begin();
    while (it != residue_ids.end()) {
        auto current = it++;
        if (rns.count(residues[*current].name()) != 0) {
            residue_ids.erase(current);
        }
    }
}

//! Remove residues which do **not** interact with a given set of other residues
//!
//! This function is designed to remove residues which do not have a desired
//! interaction with the surrounding protein environment. For example, if a user
//! is interested in small molecules that interact with a Heme group, they can
//! use this function to remove all residues that do have this interaction.
//! \param [in] frame The `frame` containing residues of interest.
//! \param [in,out] residue_ids_of_interest The residue IDs to be pruned.
//! \param [in] residue_ids_to_check The residue ids that the users wishes the 
//!  residue_ids_of_interest to interact with.
//! \param [in] distance_cutoff The distance that the residue_ids_of_interest
//!  must be within a checked residue to be included.
inline void keep_interactions(const chemfiles::Frame& frame,
                       std::set<size_t>& residue_ids_of_interest,
                       const std::set<size_t>& residue_ids_to_check,
                       double distance_cutoff = 6.0) {
    const auto& topo = frame.topology();
    const auto& residues = topo.residues();

    auto it = residue_ids_of_interest.begin();
    while (it != residue_ids_of_interest.end()) {
        auto current = it++;
        const auto& ligand_residue = residues[*current];

        for (size_t residue_to_check : residue_ids_to_check) {
            const auto& residue = residues[residue_to_check];

            for (auto prot_atom : residue) {
                for (auto lig_atom : ligand_residue) {
                    if (frame.distance(prot_atom, lig_atom) < distance_cutoff) {
                        goto found_interaction;
                    }
                }
            }
        }

        residue_ids_of_interest.erase(current);
    found_interaction:;
    }
}

//! Remove residues which **do interact** with a given set of other residues
//!
//! This function is designed to remove residues which have a undesirable
//! interaction with the surrounding protein environment. For example, if a user
//! is interested in small molecules that do not interact with water, they can
//! use this function to remove all residues that interact with water.
//! \param [in] frame The `frame` containing residues of interest.
//! \param [in,out] residue_ids_of_interest The residue IDs to be pruned.
//! \param [in] residue_ids_to_check The residue ids that the users wishes the 
//!  residue_ids_of_interest to **not** interact with.
//! \param [in] distance_cutoff The distance that the residue_ids_of_interest
//!  must be within a checked residue to be removed.
inline void remove_interactions(const chemfiles::Frame& frame,
                       std::set<size_t>& residue_ids_of_interest,
                       const std::set<size_t>& residue_ids_to_check,
                       double distance_cutoff = 6.0) {
    const auto& topo = frame.topology();
    const auto& residues = topo.residues();

    auto it = residue_ids_of_interest.begin();
    while (it != residue_ids_of_interest.end()) {
        auto current = it++;
        const auto& ligand_residue = residues[*current];

        for (size_t residue_to_check : residue_ids_to_check) {
            const auto& residue = residues[residue_to_check];

            for (auto prot_atom : residue) {
                for (auto lig_atom : ligand_residue) {
                    if (frame.distance(prot_atom, lig_atom) < distance_cutoff) {
                        residue_ids_of_interest.erase(current);
                        goto found_interaction;
                    }
                }
            }
        }
        
    found_interaction:;
    }
}

} // namespace prune

} // namespace lemon

#endif
