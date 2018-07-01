#include "benchmarker/separate.hpp"

#include <set>
#include <iostream>

using namespace benchmarker;

void benchmarker::separate_protein_and_ligand(
    const chemfiles::Frame& input,
    const chemfiles::Residue& ligand_residue,
    chemfiles::Frame& protein,
    chemfiles::Frame& ligand,
    double pocket_size
) {
    const auto& topo = input.topology();
    const auto& positions = input.positions();

    std::set<size_t> accepted_atoms;

    for (size_t prot_atom = 0; prot_atom < input.size(); ++prot_atom) {
        if (ligand_residue.contains(prot_atom)) {
            continue;
        }

        for (auto lig_atom : ligand_residue) {
            if (input.distance(prot_atom, lig_atom) < pocket_size ) {
                accepted_atoms.insert(prot_atom);
            }
        }
    }

    std::set<std::pair<size_t, size_t>> added_res;
    for (auto atom : accepted_atoms) {
        const auto& res = topo.residue_for_atom(atom);

        auto chainid = static_cast<size_t>(res->get("chainindex")->as_double());

        if (added_res.count({*(res->id()), chainid})) {
            continue;
        }

        added_res.insert({*(res->id()), chainid});

        auto res_new = chemfiles::Residue(res->name(), *(res->id()));

        for (size_t res_atom : *res) {
            protein.add_atom(input[res_atom], positions[res_atom]);
            res_new.add_atom(protein.size() - 1);
        }

        res_new.set("chainid", res->get("chainid")->as_string());
        protein.add_residue(res_new);
    }

    std::unordered_map<size_t, size_t> old_to_new;
    for (auto lig_atom : ligand_residue) {
        ligand.add_atom(input[lig_atom], positions[lig_atom]);
        old_to_new.insert({lig_atom, ligand.size() - 1});
    }

    const auto& old_bonds = topo.bonds();
    for (size_t bond_idx = 0; bond_idx < old_bonds.size(); ++bond_idx) {
        if (ligand_residue.contains(old_bonds[bond_idx][0]) && ligand_residue.contains(old_bonds[bond_idx][1])) {
            ligand.add_bond(
                old_to_new[old_bonds[bond_idx][0]],
                old_to_new[old_bonds[bond_idx][1]],
                topo.bond_orders()[bond_idx]
            );
        }
    }
}
