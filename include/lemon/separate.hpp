#ifndef LEMON_SEPARATE_HPP
#define LEMON_SEPARATE_HPP

#include <list>
#include <set>
#include <unordered_set>

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#include <chemfiles/Frame.hpp>
#include <chemfiles/Topology.hpp>
LEMON_EXTERNAL_FILE_POP

namespace lemon {

//! Separate entry frames into corresponding sub-frames
namespace separate {

//! Copy residues from one frame into a second frame
//!
//! Once a set of residues has been selected (and pruned), the residues of
//! interest should be copied into a new `frame` so that they can written to
//! disk separately. This function performs this action while taking care to
//! copy the original connectivity of the residue.
//! \param [in] input The original Frame from where the residues will be copied.
//! \param [in] accepted_residues The residue IDs for the residues to be copied.
//! \param [in,out] new_frame The frame where the residues wil be copied to.
template <typename Container>
inline void residues(const chemfiles::Frame& input,
                     const Container& accepted_residues,
                     chemfiles::Frame& new_frame) {

    const auto& residues = input.topology().residues();
    const auto& positions = input.positions();
    const auto& old_bonds = input.topology().bonds();
    const auto& bond_ord = input.topology().bond_orders();

    std::unordered_map<size_t, size_t> old_to_new;
    std::unordered_set<size_t> accepted_atoms;
    for (auto res_id : accepted_residues) {
        const auto& res = residues[res_id];

        auto res_new = chemfiles::Residue(res.name(), *(res.id()));

        for (size_t res_atom : res) {
            new_frame.add_atom(input[res_atom], positions[res_atom]);
            res_new.add_atom(new_frame.size() - 1);
            old_to_new.insert({res_atom, new_frame.size() - 1});
            accepted_atoms.insert(res_atom);
        }

        auto chainID = res.get("chainid");
        auto newID = (chainID && chainID->kind() == chemfiles::Property::STRING)?
            chainID->as_string() : "Z";

        res_new.set("chainid", newID);
        new_frame.add_residue(std::move(res_new));
    }

    for (size_t bond_idx = 0; bond_idx < old_bonds.size(); ++bond_idx) {
        if (accepted_atoms.count(old_bonds[bond_idx][0]) &&
            accepted_atoms.count(old_bonds[bond_idx][1])) {

            new_frame.add_bond(old_to_new[old_bonds[bond_idx][0]],
                               old_to_new[old_bonds[bond_idx][1]],
                               bond_ord[bond_idx]);
        }
    }
}

//! Separate a ligand and surrounding protein pocket into individual frames.
//!
//! The environment surrounding a ligand in a protein defines the *environment*
//! of that ligand. This function is meant to separate a ligand and relevent
//! *environment* into separate *frame*s s that they can written to disk.
//! \param [in] input The original Frame from where the residues will be copied.
//! \param [in] ligand_id The residue IDs for the ligand.
//! \param [in] pocket_size The radius of the ligand environment copied into
//! protein \param [in,out] protein The frame where the protein residues wil be
//! copied to. \param [in,out] ligand The frame where the ligand residue will be
//! copied to.
inline void protein_and_ligand(const chemfiles::Frame& input, size_t ligand_id,
                               double pocket_size, chemfiles::Frame& protein,
                               chemfiles::Frame& ligand) {
    const auto& topology = input.topology();
    const auto& residues = topology.residues();
    const auto& ligand_residue = residues[ligand_id];

    std::list<size_t> accepted_residues;
    for (size_t res_id = 0; res_id < residues.size(); ++res_id) {
        if (res_id == ligand_id) {
            continue;
        }

        const auto& res = residues[res_id];
        for (auto prot_atom : res) {
            for (auto lig_atom : ligand_residue) {
                if (input.distance(prot_atom, lig_atom) < pocket_size) {
                    accepted_residues.push_back(res_id);
                    goto found_interaction;
                }
            }
        }
    found_interaction:;
    }

    lemon::separate::residues(input, accepted_residues, protein);
    lemon::separate::residues(input, std::list<size_t>({ligand_id}), ligand);

    ligand.set("name", ligand_residue.name());
}

} // namespace separate

} // namespace lemon

#endif
