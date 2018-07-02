#ifndef SELECT_HPP
#define SELECT_HPP

#include <set>

#include "chemfiles/Frame.hpp"

namespace benchmarker {
    std::set<size_t> select_small_molecule(const chemfiles::Frame& input, size_t min_atoms = 10);
    std::set<size_t> select_metal_ions(const chemfiles::Frame& input);
}

#endif
