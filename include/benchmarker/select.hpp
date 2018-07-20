#ifndef SELECT_HPP
#define SELECT_HPP

#include <set>

#include "chemfiles/Frame.hpp"

namespace benchmarker
{
std::set<size_t> select_small_molecule(const chemfiles::Frame& input, size_t min_atoms = 10)
{
    const auto& residues = input.topology().residues();

    std::set<size_t> selected_residues;
    std::unordered_map<std::string, size_t> named_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size(); ++selected_residue)
    {
        const auto& residue = residues[selected_residue];

        // Quick check, filters out many ions and small molecules before more expensive checks.
        // For example, water
        if (residue.size() < min_atoms)
        {
            continue;
        }

        auto composition_type = residue.get("composition_type")->as_string();
        if (composition_type != "NON-POLYMER" && composition_type != "OTHER" && composition_type != "PEPTIDE-LIKE")
        {
            continue;
        }

        size_t num_heavy_atoms = std::count_if(
            std::begin(residue),
            std::end(residue),
            [&input](size_t index) {
                return *(input[index].atomic_number()) != 1;
            });

        if (num_heavy_atoms < min_atoms)
        {
            continue;
        }

        selected_residues.insert(selected_residue);
    }

    return selected_residues;
}

std::set<size_t> select_metal_ions(const chemfiles::Frame& input)
{
    const auto& residues = input.topology().residues();

    std::set<size_t> selected_residues;
    std::unordered_map<std::string, size_t> named_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size(); ++selected_residue)
    {
        const auto& residue = residues[selected_residue];

        if (residue.size() == 1 && input[*residue.begin()].charge() > 0.0)
        {
            selected_residues.insert(selected_residue);
        }
    }

    return selected_residues;
}
}

#endif
