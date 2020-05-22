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
//! \param [in] entry The original Frame from where the residues will be copied.
//! \param [in] accepted_residues The residue IDs for the residues to be copied.
//! \param [in,out] new_frame The frame where the residues wil be copied to.
template <typename Container>
inline void residues(const chemfiles::Frame& entry,
                     const Container& accepted_residues,
                     chemfiles::Frame& new_frame,
                     const std::string& use_altloc = "A") {

    const auto& residues = entry.topology().residues();
    const auto& positions = entry.positions();
    const auto& old_bonds = entry.topology().bonds();
    const auto& bond_ord = entry.topology().bond_orders();

    std::unordered_map<size_t, size_t> old_to_new;
    std::unordered_set<size_t> accepted_atoms;
    for (auto res_id : accepted_residues) {
        const auto& res = residues[res_id];

        auto res_new = chemfiles::Residue(res.name(), *(res.id()));

        for (size_t res_atom : res) {

            auto altloc = entry[res_atom].get<chemfiles::Property::STRING>("altloc").value_or(" ");
            if (altloc != " " && altloc != use_altloc) {
                continue;
            }

            new_frame.add_atom(entry[res_atom], positions[res_atom]);
            res_new.add_atom(new_frame.size() - 1);
            old_to_new.insert({res_atom, new_frame.size() - 1});
            accepted_atoms.insert(res_atom);
        }

        for (const auto& prop : res.properties()) {
            if (prop.first == "chainid") {
                res_new.set(prop.first, std::string(1, prop.second.as_string()[0]));
                continue;
            }
            res_new.set(prop.first, prop.second);
        }

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
//! \param [in] entry The original Frame from where the residues will be copied.
//! \param [in] ligand_id The residue ID for the ligand.
//! \param [in] pocket_size The radius of the ligand environment copied into
//! protein
//! \param [in,out] protein The frame where the protein residues wil be
//! copied to.
//! \param [in,out] ligand The frame where the ligand residue will be
//! copied to.
inline void protein_and_ligand(const chemfiles::Frame& entry, size_t ligand_id,
                               double pocket_size, chemfiles::Frame& protein,
                               chemfiles::Frame& ligand,
                               const std::string& altloc = "A"
                               ) {
    const auto& topology = entry.topology();
    const auto& residues = topology.residues();
    const auto& ligand_residue = residues[ligand_id];

    std::list<size_t> accepted_residues;
    for (size_t res_id = 0; res_id < residues.size(); ++res_id) {
        if (res_id == ligand_id) {
            continue;
        }

        const auto& res = residues[res_id];

        // Do some cleaning up here
        if (res.name() == "UNX" || res.name() == "UNL") {
            continue;
        }

        for (auto prot_atom : res) {
            for (auto lig_atom : ligand_residue) {
                if (pocket_size <= 0 ||
                    entry.distance(prot_atom, lig_atom) < pocket_size) {
                    accepted_residues.push_back(res_id);
                    goto found_interaction;
                }
            }
        }
    found_interaction:;
    }

    lemon::separate::residues(entry, accepted_residues, protein, altloc);
    lemon::separate::residues(entry, std::list<size_t>({ligand_id}), ligand, altloc);

    ligand.set("name", ligand_residue.name());
}

//! Separate ligands and surrounding protein pocket into individual frames.
//!
//! The environment surrounding ligands in a protein defines the *environment*
//! of those ligands. This function is meant to separate ligands and relevent
//! *environment* into separate *frame*s s that they can written to disk.
//! \param [in] entry The original Frame from where the residues will be copied.
//! \param [in] ligand_ids The residue IDs for the ligand.
//! \param [in] pocket_size The radius of the ligand environment copied into
//! protein
//! \param [in,out] protein The frame where the protein residues wil be
//! copied to.
//! \param [in,out] ligand The frame where the ligand residue will be
//! copied to.
template <typename Container>
inline void protein_and_ligands(const chemfiles::Frame& entry,
                                const Container& ligand_ids,
                                double pocket_size, chemfiles::Frame& protein,
                                chemfiles::Frame& ligand,
                                const std::string& altloc = "A") {
    const auto& topology = entry.topology();
    const auto& residues = topology.residues();

    std::list<size_t> accepted_residues;
    for (size_t res_id = 0; res_id < residues.size(); ++res_id) {
        const auto& res = residues[res_id];

        // Do some cleaning up here
        if (res.name() == "UNX" || res.name() == "UNL") {
            continue;
        }

        for (auto lig_res : ligand_ids) {
            if (lig_res == res_id) {
                goto found_interaction;
            }
        }

        for (auto lig_res : ligand_ids) {
            for (auto prot_atom : res) {
                for (auto lig_atom : residues[lig_res]) {
                    if (pocket_size <= 0 ||
                        entry.distance(prot_atom, lig_atom) < pocket_size) {
                        accepted_residues.push_back(res_id);
                        goto found_interaction;
                    }
                }
            }
        }
    found_interaction:;
    }

    lemon::separate::residues(entry, accepted_residues, protein);
    lemon::separate::residues(entry, ligand_ids, ligand, altloc);
}

} // namespace separate

} // namespace lemon

#endif
