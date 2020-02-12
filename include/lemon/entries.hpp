#ifndef LEMON_ENTRIES_HPP
#define LEMON_ENTRIES_HPP

#include <array>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <cctype>
#include <unordered_set>
#include <vector>
#include <algorithm>

#include "lemon/residue_name.hpp"

namespace lemon {

using Entries = std::unordered_set<std::string>;

inline std::string::value_type toupper(std::string::value_type ch) {
    return static_cast<std::string::value_type>(std::toupper(ch));
}

inline std::string::value_type tolower(std::string::value_type ch) {
    return static_cast<std::string::value_type>(std::tolower(ch));
}

inline Entries read_entry_file(std::istream& input) {
    Entries result;
    std::string temp;
    while (std::getline(input, temp)) {
        std::string a = temp.substr(0, 4);
        std::transform(a.begin(), a.end(), a.begin(), toupper);
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
    std::string temp;
    std::string item;
    while (std::getline(input, temp)) {
        std::string a = temp.substr(0, 4);

        std::transform(a.begin(), a.end(), a.begin(), toupper);

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
