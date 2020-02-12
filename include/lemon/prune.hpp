#ifndef LEMON_PRUNE_HPP
#define LEMON_PRUNE_HPP

#include <algorithm>
#include <list>

#include "lemon/residue_name.hpp"

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#include <chemfiles.hpp>
LEMON_EXTERNAL_FILE_POP

namespace lemon {

//! Prune selected residues by removing them based on a criterion
namespace prune {

//! The default distance used for pruning
auto constexpr DEFAULT_DISTANCE = 6.0;

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
template <typename Container>
inline void identical_residues(const chemfiles::Frame& frame,
                               Container& residue_ids) {
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
                auto bio_current = res_current.get("assembly");
                auto bio_check = res_check.get("assembly");

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
//! as sugars and fatty acids used to induce crystallization. As a result, some
//! users may remove these cofactors as they may match other criteria (such as
//! being a small molecule) set by the user.
//! \param [in] frame The `frame` containing residues of interest.
//! \param [in,out] residue_ids The residue IDs to be pruned.
//! \param [in] rns The residue names that one wishes to remove from
//! residue_ids.
template <typename Container>
inline void cofactors(const chemfiles::Frame& frame, Container& residue_ids,
                      const ResidueNameSet& rns) {
    auto& residues = frame.topology().residues();

    // In C++-17 there's a better version called std::erase_if
    residue_ids.erase(
        std::remove_if(residue_ids.begin(), residue_ids.end(),
                       [&residues, &rns](size_t current) {
                           return rns.count(residues[current].name()) != 0;
                       }),
        residue_ids.end());
}

template <typename Container1, typename Container2 = Container1>
inline void interactions(const chemfiles::Frame& frame, Container1& residue_ids,
                         const Container2& interaction_ids,
                         double distance_cutoff = DEFAULT_DISTANCE,
                         bool keep = true) {
    const auto& residues = frame.topology().residues();

    residue_ids.erase(
        std::remove_if(residue_ids.begin(), residue_ids.end(),
                       [&frame, &residues, &interaction_ids, distance_cutoff,
                        keep](size_t current) {
                           const auto& ligand_residue = residues[current];

                           for (auto residue_to_check : interaction_ids) {
                               const auto& residue = residues[residue_to_check];

                               for (auto prot_atom : residue) {
                                   for (auto lig_atom : ligand_residue) {
                                       if (frame.distance(prot_atom, lig_atom) <
                                           distance_cutoff) {
                                           return !keep;
                                       }
                                   }
                               }
                           }

                           return keep;
                       }),
        residue_ids.end());
}

//! Remove residues which do **not** interact with a given set of other residues
//!
//! This function is designed to remove residues which do not have a desired
//! interaction with the surrounding protein environment. For example, if a user
//! is interested in small molecules that interact with a Heme group, they can
//! use this function to remove all residues that do have this interaction.
//! \param [in] frame The `frame` containing residues of interest.
//! \param [in,out] residue_ids The residue IDs to be pruned.
//! \param [in] interaction_ids The residue ids that the users wishes the
//!  residue_ids to interact with.
//! \param [in] distance_cutoff The distance that the residue_ids
//!  must be within a checked residue to be included.
template <typename Container1, typename Container2 = Container1>
inline void keep_interactions(const chemfiles::Frame& frame,
                              Container1& residue_ids,
                              const Container2& interaction_ids,
                              double distance_cutoff = DEFAULT_DISTANCE) {
    interactions(frame, residue_ids, interaction_ids, distance_cutoff, true);
}

//! Remove residues which **do interact** with a given set of other residues
//!
//! This function is designed to remove residues which have a undesirable
//! interaction with the surrounding protein environment. For example, if a user
//! is interested in small molecules that do not interact with water, they can
//! use this function to remove all residues that interact with water.
//! \param [in] frame The `frame` containing residues of interest.
//! \param [in,out] residue_ids The residue IDs to be pruned.
//! \param [in] interaction_ids The residue ids that the users wishes the
//!  residue_ids to **not** interact with.
//! \param [in] distance_cutoff The distance that the residue_ids
//!  must be within a checked residue to be removed.
template <typename Container1, typename Container2 = Container1>
inline void remove_interactions(const chemfiles::Frame& frame,
                                Container1& residue_ids,
                                const Container2& interaction_ids,
                                double distance_cutoff = DEFAULT_DISTANCE) {
    interactions(frame, residue_ids, interaction_ids, distance_cutoff, false);
}

} // namespace prune

} // namespace lemon

#endif
