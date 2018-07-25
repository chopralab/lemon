#ifndef ENTRIES_HPP
#define ENTRIES_HPP

#include <array>
#include <istream>
#include <string>
#include <vector>
#include <fstream>

#include "lemon/residue_name.hpp"
#include "chemfiles/utils.hpp"

namespace lemon {

typedef std::vector<std::array<char, 4>> PDBIDVec;

void read_entry_file(std::istream& input,
                     PDBIDVec& result) {
    std::string temp;
    while (std::getline(input, temp)) {
        std::array<char, 4> a;
        a[0] = temp[0];
        a[1] = temp[1];
        a[2] = temp[2];
        a[3] = temp[3];
        result.emplace_back(a);
    }
}
void read_entry_file(const std::string& input,
                     PDBIDVec& result,
                     size_t number_of_entries = 143000) {
    std::ifstream input_file(input);
    std::string header;
    std::string junk;
    std::getline(input_file, header);
    std::getline(input_file, junk);

    result.reserve(number_of_entries);
    lemon::read_entry_file(input_file, result);
}

void read_entry_file(std::istream& input,
                     PDBIDVec& result,
                     std::unordered_map<std::string, ResidueNameSet>& rnm) {
    std::string temp;
    while (std::getline(input, temp)) {
        std::array<char, 4> a;
        a[0] = temp[0];
        a[1] = temp[1];
        a[2] = temp[2];
        a[3] = temp[3];

        const std::string pdbid = std::string(a.data(), 4);
        auto residue_names = chemfiles::split(temp, '\t');
        for (size_t i = 1; i < residue_names.size(); i += 2) {
            rnm[pdbid].insert(residue_names[i]);
        }

        result.emplace_back(a);
    }
}
}

#endif
