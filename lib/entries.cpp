#include "benchmarker/entries.hpp"

#include <istream>
#include <iostream>

void benchmarker::read_entry_file(std::istream& input, std::vector<std::array<char, 4>>& result) {
    std::string temp;
    while ( std::getline(input,temp) ) {
        std::array<char, 4> a;
        a[0] = temp[0];
        a[1] = temp[1];
        a[2] = temp[2];
        a[3] = temp[3];
        result.emplace_back(a);
    }
}
