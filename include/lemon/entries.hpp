#ifndef LEMON_ENTRIES_HPP
#define LEMON_ENTRIES_HPP

#include <array>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "lemon/residue_name.hpp"

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

    if (input == "") {
        return Entries();
    }

    std::ifstream input_file(input);
    std::string header;
    std::string junk;
    std::getline(input_file, header);
    std::getline(input_file, junk);

    return lemon::read_entry_file(input_file);
}

inline void
read_entry_file(std::istream& input, Entries& result,
                std::unordered_map<std::string, ResidueNameSet>& rnm) {
    std::string temp, item;
    while (std::getline(input, temp)) {
        std::string a = temp.substr(0, 4);

        std::stringstream ss(temp);
        std::vector<std::string> residue_names;
        size_t count = 0;
        while (std::getline(ss, item, '\t')) {
            if (++count % 2 == 0) { // check if even
                rnm[a].insert(item);
            }
        }

        result.insert(a);
    }
}
} // namespace lemon

#endif
