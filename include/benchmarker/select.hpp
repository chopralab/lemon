#ifndef SELECT_HPP
#define SELECT_HPP

#include "chemfiles/Frame.hpp"

namespace benchmarker {
    const chemfiles::Residue& select_small_molecule(const chemfiles::Frame& input, size_t min_atoms = 10);
}

#endif
