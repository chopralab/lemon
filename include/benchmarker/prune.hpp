
#ifndef PRUNE_HPP
#define PRUNE_HPP

#include <chemfiles.hpp>
#include "benchmarker/residue_name.hpp"

namespace benchmarker {
const ResidueNameSet common_cofactors{
    {"FAD"}, {"FMN"}, {"NAD"}, {"NAP"}, {"CLA"}, {"HEM"}, {"HEA"},
    {"HEB"}, {"HEC"}, {"ADP"}, {"ATP"}, {"GDP"}, {"GTP"}, {"UNL"},
};

void remove_identical_residues(const chemfiles::Frame& file,
                               std::set<size_t>& residue_ids) {
    auto& residues = file.topology().residues();

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

void remove_common_cofactors(const chemfiles::Frame& file,
                             std::set<size_t>& residue_ids) {
    auto& residues = file.topology().residues();

    auto it = residue_ids.begin();
    while (it != residue_ids.end()) {
        auto current = it++;
        if (common_cofactors.count(residues[*current].name()) != 0) {
            residue_ids.erase(current);
        }
    }
}

bool find_nucleic_acid_interactions(const chemfiles::Frame& input,
                                    std::set<size_t>& residue_ids,
                                    double distance_cutoff = 6.0) {
    const auto& topo = input.topology();
    const auto& positions = input.positions();
    const auto& residues = topo.residues();

    auto it = residue_ids.begin();
    while (it != residue_ids.end()) {
        auto current = it++;
        const auto& ligand_residue = residues[*current];

        bool has_interaction = false;
        for (size_t prot_atom = 0; prot_atom < input.size(); ++prot_atom) {
            if (ligand_residue.contains(prot_atom)) {
                continue;
            }

            for (auto lig_atom : ligand_residue) {
                if (input.distance(prot_atom, lig_atom) < distance_cutoff) {
                    auto comp_type = topo.residue_for_atom(prot_atom)
                                         ->get("composition_type")
                                         ->as_string();
                    if (comp_type.find("DNA") != std::string::npos ||
                        comp_type.find("RNA") != std::string::npos) {
                            has_interaction = true;
                            break;
                    }
                }
            }

            if (has_interaction) {
                break;
            }
        }

        if (!has_interaction) {
            residue_ids.erase(current);
        }
    }
}
}

#endif
