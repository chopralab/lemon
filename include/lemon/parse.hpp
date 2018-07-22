#ifndef PARSE_HPP
#define PARSE_HPP

#include <set>
#include <sstream>

#include <chemfiles.hpp>
#include <chemfiles/utils.hpp>
#include "lemon/entries.hpp"
#include "lemon/residue_name.hpp"

namespace lemon {
void count_residues(const chemfiles::Frame& file,
                    ResidueNameCount& resn_count) {
    auto& residues = file.topology().residues();

    for (auto& residue : residues) {
        auto resn_count_iterator = resn_count.find(residue.name());

        if (resn_count_iterator == resn_count.end()) {
            resn_count[residue.name()] = 1;
            continue;
        }

        ++resn_count_iterator->second;
    }
}

size_t count_bioassemblies(const chemfiles::Frame& file) {
    auto& residues = file.topology().residues();
    std::set<std::string> assembies;

    for (auto& residue : residues) {
        if (residue.get("assembly")) {
            assembies.insert(residue.get("assembly")->as_string());
        }
    }

    return assembies.size();
}

void print_residue_name_counts(std::ostream& os, const std::string& pdbid,
                               const chemfiles::Frame& complex,
                               const std::set<size_t>& res_ids) {
    if (res_ids.empty()) {
        return;
    }

    std::unordered_map<std::string, size_t> residue_counts;
    const auto& residues = complex.topology().residues();
    for (auto res_id : res_ids) {
        auto iter = residue_counts.find(residues[res_id].name());
        if (iter == residue_counts.end()) {
            residue_counts[residues[res_id].name()] = 1;
            continue;
        }
        ++iter->second;
    }

    std::stringstream ss;
    ss << pdbid;
    for (const auto iter : residue_counts) {
        ss << " " << iter.first << " " << iter.second;
    }
    ss << "\n";

    os << ss.str();
}
}

#endif
