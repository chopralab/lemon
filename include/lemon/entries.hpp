#ifndef LEMON_ENTRIES_HPP
#define LEMON_ENTRIES_HPP

#include <array>
#include <istream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_set>

#include "lemon/residue_name.hpp"
#include "chemfiles/utils.hpp"

namespace lemon {

typedef std::unordered_set<std::string> Entries;

inline Entries read_entry_file(std::istream& input) {
    Entries result;
    std::string temp;
    while (std::getline(input, temp)) {
        std::string a = temp.substr(0, 4);
        result.insert(a);
    }
    return result;
}
inline Entries read_entry_file(const std::string& input) {
    std::ifstream input_file(input);
    std::string header;
    std::string junk;
    std::getline(input_file, header);
    std::getline(input_file, junk);

    return lemon::read_entry_file(input_file);
}

inline void read_entry_file(std::istream& input,
                     Entries& result,
                     std::unordered_map<std::string, ResidueNameSet>& rnm) {
    std::string temp;
    while (std::getline(input, temp)) {
        std::string a = temp.substr(0, 4);

        const std::string pdbid = std::string(a.data(), 4);
        auto residue_names = chemfiles::split(temp, '\t');
        for (size_t i = 1; i < residue_names.size(); i += 2) {
            rnm[pdbid].insert(residue_names[i]);
        }

        result.insert(a);
    }
}
}

#endif
