
#ifndef PARSE_HPP
#define PARSE_HPP

#include <chemfiles.hpp>
#include "benchmarker/residue_name.hpp"

namespace benchmarker {
void retreive_residue_counts(const chemfiles::Frame& file, ResidueNameCount& resn_count) {
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
}

#endif
