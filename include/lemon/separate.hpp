#ifndef SEPARATE_HPP
#define SEPARATE_HPP

#include <set>

#include "chemfiles/Topology.hpp"
#include "chemfiles/Frame.hpp"

namespace lemon {
inline void separate_protein_and_ligand(const chemfiles::Frame& input,
                                 size_t ligand_id,
                                 chemfiles::Frame& protein,
                                 chemfiles::Frame& ligand, double pocket_size) {
    const auto& topo = input.topology();
    const auto& positions = input.positions();
    const auto& residues = topo.residues();
    const auto& ligand_residue = residues[ligand_id];

    std::set<size_t> accepted_residues;
    for (size_t res_id = 0; res_id < residues.size(); ++res_id) {
        if (res_id == ligand_id) {
            continue;
        }

        const auto& res = residues[res_id];
        for (auto prot_atom : res) {
            for (auto lig_atom : ligand_residue) {
                if (input.distance(prot_atom, lig_atom) < pocket_size) {
                    accepted_residues.insert(res_id);
                    goto found_interaction;
                }
            }
        }
        found_interaction:;
    }

    for (auto res_id : accepted_residues) {
        const auto& res = residues[res_id];

        auto res_new = chemfiles::Residue(res.name(), *(res.id()));

        for (size_t res_atom : res) {
            protein.add_atom(input[res_atom], positions[res_atom]);
            res_new.add_atom(protein.size() - 1);
        }

        res_new.set("chainid", res.get("chainid")->as_string());
        protein.add_residue(res_new);
    }

    std::unordered_map<size_t, size_t> old_to_new;
    for (auto lig_atom : ligand_residue) {
        ligand.add_atom(input[lig_atom], positions[lig_atom]);
        old_to_new.insert({lig_atom, ligand.size() - 1});
    }

    const auto& old_bonds = topo.bonds();
    for (size_t bond_idx = 0; bond_idx < old_bonds.size(); ++bond_idx) {
        if (ligand_residue.contains(old_bonds[bond_idx][0]) &&
            ligand_residue.contains(old_bonds[bond_idx][1])) {
            ligand.add_bond(old_to_new[old_bonds[bond_idx][0]],
                            old_to_new[old_bonds[bond_idx][1]],
                            topo.bond_orders()[bond_idx]);
        }
    }

    ligand.set("name", ligand_residue.name());
}
}

#endif
