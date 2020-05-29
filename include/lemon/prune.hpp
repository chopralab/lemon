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
inline Container identical_residues(const chemfiles::Frame& frame,
                                    Container& residue_ids) {
    auto& residues = frame.topology().residues();

    if (residue_ids.empty()) {
        return residue_ids;
    }

    auto assembly = residues[*residue_ids.begin()].get("assembly");

    residue_ids.erase(
        std::remove_if(residue_ids.begin(), residue_ids.end(),
            [&residues, assembly](size_t current_id) {
                auto& residue = residues[current_id];
                auto prop = residue.get("assembly");

                return prop != assembly;
            }
        ),
    residue_ids.end());

    return residue_ids;
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
inline Container cofactors(const chemfiles::Frame& frame, Container& residue_ids,
                           const ResidueNameSet& rns) {
    auto& residues = frame.topology().residues();

    // In C++-17 there's a better version called std::erase_if
    residue_ids.erase(
        std::remove_if(residue_ids.begin(), residue_ids.end(),
                       [&residues, &rns](size_t current) {
                           return rns.count(residues[current].name()) != 0;
                       }),
        residue_ids.end()
    );

    return residue_ids;
}

template <typename Container1, typename Container2 = Container1>
inline Container1 interactions(const chemfiles::Frame& frame,
                              Container1& residue_ids,
                              const Container2& interaction_ids,
                              double distance_cutoff = DEFAULT_DISTANCE,
                              bool keep = true) {
    const auto& residues = frame.topology().residues();

    residue_ids.erase(
        std::remove_if(residue_ids.begin(), residue_ids.end(),
            [&](size_t current) {
                const auto& ligand_residue = residues[current];
            
                for (auto residue_to_check : interaction_ids) {
                    const auto& residue = residues[residue_to_check];
            
                    for (auto prot_atom : residue) {
                        for (auto lig_atom : ligand_residue) {
                            if (distance_cutoff <= 0 || 
                                frame.distance(prot_atom, lig_atom) <
                                distance_cutoff) {
                                return !keep;
                            }
                        }
                    }
                }
            
                return keep;
            }),
        residue_ids.end()
    );

    return residue_ids;
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
inline Container1 keep_interactions(const chemfiles::Frame& frame,
                                    Container1& residue_ids,
                                    const Container2& interaction_ids,
                                    double distance_cutoff = DEFAULT_DISTANCE) {
    return interactions(frame, residue_ids, interaction_ids, distance_cutoff, true);
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
inline Container1 remove_interactions(const chemfiles::Frame& frame,
                                      Container1& residue_ids,
                                      const Container2& interaction_ids,
                                      double distance_cutoff = DEFAULT_DISTANCE) {
    return interactions(frame, residue_ids, interaction_ids, distance_cutoff, false);
}

//! Turns `residue_ids` in to intersection between it and `intersection_ids`
//!
//! This function is designed to keep residues which have a desirable
//! intersection with another set of residue ids.
//! \param [in,out] residue_ids The residue IDs to be pruned.
//! \param [in] intersection_ids The residue ids that the users wishes the
//!  residue_ids to also be in.
template <typename Container1, typename Container2 = Container1>
inline Container1 intersection(Container1& residue_ids,
                               const Container2& intersection_ids) {

    // ugly O(n*m), but gets the job done. Ideally one should use sets though.
    residue_ids.erase(
        std::remove_if(residue_ids.begin(), residue_ids.end(),
            [&intersection_ids](size_t current){
                for (auto& i : intersection_ids){
                    if (i == current) {
                        return false;
                    }
                }
                return true;
        }), residue_ids.end()
    );

    return residue_ids;
}

//! Keeps residues with a given property
//!
//! This function is designed to keep residues which have a desirable property
//! \param [in] frame The `frame` containing residues of interest.
//! \param [in,out] residue_ids The residue IDs to be pruned.
//! \param [in] property_name The name of the property to keep
//! \param [in] property The property that the residues must have to be kept
template <typename Container>
inline Container has_property(const chemfiles::Frame& frame,
                               Container& residue_ids,
                               const std::string& property_name,
                               const chemfiles::Property& property) {

    residue_ids.erase(
        std::remove_if(residue_ids.begin(), residue_ids.end(),
            [&frame, &property_name, &property](size_t current){
                auto& residue = frame.topology().residues()[current];
                auto prop = residue.get(property_name);
                if (!prop) {
                    return true;
                }

                return *prop != property;
        }), residue_ids.end()
    );

    return residue_ids;
}

} // namespace prune

} // namespace lemon

#endif
