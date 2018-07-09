
#ifndef PRUNE_HPP
#define PRUNE_HPP

#include <chemfiles.hpp>
#include "benchmarker/residue_name.hpp"

namespace benchmarker {
    void remove_identical_residues(const chemfiles::Frame& file, std::set<size_t>& residue_ids);
    void remove_common_cofactors(const chemfiles::Frame& file, std::set<size_t>& residue_ids);
}

#endif
