#ifndef ENTRIES_HPP
#define ENTRIES_HPP

#include <string>
#include <vector>
#include <array>

namespace benchmarker {
    void read_entry_file(std::istream& file, std::vector<std::array<char, 4>>& result);
}

#endif
