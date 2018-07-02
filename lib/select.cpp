#include "benchmarker/select.hpp"

#include <set>

using namespace benchmarker;

const chemfiles::Residue& benchmarker::select_small_molecule(const chemfiles::Frame& input, size_t min_atoms) {
    const auto& residues = input.topology().residues();

    std::set<size_t> selected_residues;
    std::unordered_map<std::string, size_t> named_residues;

    for (size_t selected_residue = 0; selected_residue < residues.size(); ++selected_residue) {
        const auto& residue = residues[selected_residue];

        auto composition_type = residue.get("composition_type")->as_string();
        if ( composition_type != "NON-POLYMER" && composition_type != "OTHER" && composition_type != "PEPTIDE-LIKE") {
            continue;
        }

        if (residue.size() < min_atoms) {
            continue;
        }

        selected_residues.insert(selected_residue);

        auto iter = named_residues.find(residue.name());
        if (iter == named_residues.end()) {
            named_residues[residue.name()] = 0;
            continue;
        }

        iter->second++;
    }

    if (selected_residues.size() == 0) {
        throw std::range_error("No suitable residue found");
    } else if (selected_residues.size() == 1) {
        return residues[ *selected_residues.begin() ];
    }

    throw std::range_error("Multiple suitable residues found");
}
