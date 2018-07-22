#ifndef ENTRIES_HPP
#define ENTRIES_HPP

#include <array>
#include <istream>
#include <string>
#include <vector>
#include <fstream>

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
                     size_t number_of_entries = 142000) {
    std::ifstream input_file(input);
    std::string header;
    std::string junk;
    std::getline(input_file, header);
    std::getline(input_file, junk);

    result.reserve(number_of_entries);
    lemon::read_entry_file(input_file, result);
}
}

#endif
