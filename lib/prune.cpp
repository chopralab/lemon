#include "benchmarker/prune.hpp"

#include <set>

namespace benchmarker {

const ResidueNameSet common_cofactors {
    {"FAD"}, {"FMN"},
    {"NAD"}, {"NAP"},
    {"CLA"},
    {"HEM"}, {"HEA"}, {"HEB"}, {"HEC"},
    {"ADP"}, {"ATP"},
    {"GDP"}, {"GTP"}
};

void remove_identical_residues(const chemfiles::Frame& file, std::set<size_t>& residue_ids) {
    auto& residues = file.topology().residues();

    auto it = residue_ids.begin();
    while(it != residue_ids.end()) {
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

                // We use bioassemblies because they are standardized by the protein databank
                auto bio_current = res_current.get("assembly")? res_current.get("assembly")->as_string() : "";
                auto bio_check = res_check.get("assembly")? res_check.get("assembly")->as_string() : "";

                if (bio_current != bio_check) {
                    residue_ids.erase(check_current);
                }
            }
        }

        ++it;
    }
}

void remove_common_cofactors(const chemfiles::Frame& file, std::set<size_t>& residue_ids) {

    auto& residues = file.topology().residues();

    auto it = residue_ids.begin();
    while(it != residue_ids.end()) {
        auto current = it++;
        if (common_cofactors.count(residues[*current].name()) != 0) {
            residue_ids.erase(current);
        }
    }
}

}
