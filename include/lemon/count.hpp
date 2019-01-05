#ifndef LEMON_PARSE_HPP
#define LEMON_PARSE_HPP

#include <set>
#include <list>
#include <sstream>

#include <chemfiles.hpp>
#include <chemfiles/utils.hpp>
#include "lemon/entries.hpp"
#include "lemon/residue_name.hpp"

namespace lemon {

//! \brief Count various biological features
namespace count {

//! Append all residue counts in a `Frame` to a `ResidueNameCount`
//!
//! Use this function to count the all residues in `frame` using the residue
//! name.  If the residue is previously stored in the supplied `resn_count`
//! then the count of the residue is increased.
//! \param [in] frame The frame containing residues of interest
//! \param [in,out] resn_count A map of residue names to their respective count
inline void residues(const chemfiles::Frame& frame,
                     ResidueNameCount& resn_count) {
    auto& residues = frame.topology().residues();

    for (auto& residue : residues) {
        auto resn_count_iterator = resn_count.find(residue.name());

        if (resn_count_iterator == resn_count.end()) {
            resn_count[residue.name()] = 1;
            continue;
        }

        ++resn_count_iterator->second;
    }
}

//! Append selected residue counts in a `Frame` to a `ResidueNameCount`
//!
//! Use this function to count the residues in `frame` using the residue
//! name with the residue ids in `resids`. If the residue is previously stored
//! in the supplied `resn_count` then the count of the residue is increased.
//! \param [in] frame The frame containing residues of interest
//! \param [in] resids Residue ids to consider
//! \param [in,out] resn_count A map of residue names to their respective count
template<typename Container>
inline void residues(const chemfiles::Frame& frame,
                     const Container& resids,
                     ResidueNameCount& resn_count) {
    auto& residues = frame.topology().residues();

    for (auto& resid : resids) {
        auto resn_count_iterator = resn_count.find(residues[resid].name());

        if (resn_count_iterator == resn_count.end()) {
            resn_count[residues[resid].name()] = 1;
            continue;
        }

        ++resn_count_iterator->second;
    }
}

//! Obtain the number of alternative location types in a `Frame`.
//!
//! \param [in] frame The frame containing atoms of interest.
//! \return the number of unique alternative location names
inline size_t altloc(const chemfiles::Frame& frame) {
    std::set<char> alt_locs;
    for (const auto& atom : frame) {
        const auto& altloc = atom.get("altloc");
        if (altloc) {
            alt_locs.insert(altloc->as_string()[0]);
        }
    }

    return alt_locs.size();
}

//! Obtain the number of bioassemblies in a `Frame`.
//!
//! \param [in] frame The frame containing assemblies of interest.
//! \return the number of unique bioassemblies location names
inline size_t bioassemblies(const chemfiles::Frame& frame) {
    auto& residues = frame.topology().residues();
    std::set<std::string> assembies;

    for (auto& residue : residues) {
        if (residue.get("assembly")) {
            assembies.insert(residue.get("assembly")->as_string());
        }
    }

    return assembies.size();
}

//! Print select residue names and their respective counts
//!
//! \param [in] pdbid The current PDB id of the complex to be printed
//! \param [in] complex The complex containing the residues of interest
//! \param [in] res_ids The set of residue ids for printing
//! \return a string containing the formatted output.
template<typename Container>
inline std::string print_residue_name_counts(const std::string& pdbid,
                                             const chemfiles::Frame& complex,
                                             const Container& res_ids) {
    if (res_ids.empty()) {
        return "";
    }

    ResidueNameCount rnc;
    lemon::count::residues(complex, res_ids, rnc);

    std::stringstream ss;
    ss << pdbid << rnc << "\n";
    return ss.str();
}

} // namespace count

} // namespace lemon

#endif
