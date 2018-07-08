
#ifndef PARSE_HPP
#define PARSE_HPP

#include <chemfiles.hpp>
#include "benchmarker/residue_name.hpp"

namespace benchmarker {
    void retreive_residue_counts(const chemfiles::Frame&, ResidueNameCount& resn_count);
    size_t count_bioassemblies(const chemfiles::Frame&);
}

#endif
