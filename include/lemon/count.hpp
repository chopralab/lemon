#ifndef PARSE_HPP
#define PARSE_HPP

#include <set>
#include <sstream>

#include <chemfiles.hpp>
#include <chemfiles/utils.hpp>
#include "lemon/entries.hpp"
#include "lemon/residue_name.hpp"

namespace lemon {
inline void count_residues(const chemfiles::Frame& file,
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

inline void count_residues(const chemfiles::Frame& file,
                    const std::set<size_t>& resids,
                    ResidueNameCount& resn_count) {
    auto& residues = file.topology().residues();

    for (auto& resid : resids) {
        auto resn_count_iterator = resn_count.find(residues[resid].name());

        if (resn_count_iterator == resn_count.end()) {
            resn_count[residues[resid].name()] = 1;
            continue;
        }

        ++resn_count_iterator->second;
    }
}

inline size_t count_altloc(const chemfiles::Frame& files) {
    std::set<char> alt_locs;
    for (const auto& atom : files) {
        const auto& altloc = atom.get("altloc");
        if (altloc) {
            alt_locs.insert(altloc->as_string()[0]);
        }
    }

    return alt_locs.size();
}

inline size_t count_bioassemblies(const chemfiles::Frame& file) {
    auto& residues = file.topology().residues();
    std::set<std::string> assembies;

    for (auto& residue : residues) {
        if (residue.get("assembly")) {
            assembies.insert(residue.get("assembly")->as_string());
        }
    }

    return assembies.size();
}

inline void print_residue_name_counts(std::ostream& os, const std::string& pdbid,
                               const chemfiles::Frame& complex,
                               const std::set<size_t>& res_ids) {
    if (res_ids.empty()) {
        return;
    }

    ResidueNameCount rnc;
    count_residues(complex, res_ids, rnc);

    std::stringstream ss;
    ss << pdbid << rnc << "\n";
    os << ss.str();
}
}

#endif
