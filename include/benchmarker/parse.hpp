
#ifndef MMTF_PARSE_HPP
#define MMTF_PARSE_HPP

#include <chemfiles.hpp>
#include <unordered_map>

namespace benchmarker{

    typedef std::unordered_map<std::string, std::size_t> ResnCount;

    void retreive_residue_counts(const chemfiles::Frame&, ResnCount& resn_count);
}

#endif
