#ifndef LEMON_PARSE_HPP
#define LEMON_PARSE_HPP

#include <list>
#include <set>
#include <sstream>
#include <unordered_set>

#include "lemon/entries.hpp"
#include "lemon/residue_name.hpp"
#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#include <chemfiles.hpp>
LEMON_EXTERNAL_FILE_POP

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
template <typename Container>
inline void residues(const chemfiles::Frame& frame, const Container& resids,
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

//! Obtain the number of times a given property occurs in a `Frame`.
//!
//! \param [in] frame The frame containing atoms of interest.
//! \param [in] name The name of the property
//! \return the number of unique alternative location names
inline size_t atom_property(const chemfiles::Frame& frame,
                            const std::string& name) {
    std::unordered_set<std::string> string_props;
    for (const auto& atom : frame) {
        const auto& prop = atom.get(name);
        if (!prop) {
            continue;
        }
        switch (prop->kind()) {
        case chemfiles::Property::STRING:
            string_props.insert(prop->as_string());
            break;
        default:
            break; // Not a type we support
        }
    }

    return string_props.size();
}

//! Obtain the number of bioassemblies in a `Frame`.
//!
//! \param [in] frame The frame containing residues of interest.
//! \param [in] name The name of the property
//! \return the number of unique bioassemblies location names
inline size_t residue_property(const chemfiles::Frame& frame,
                               const std::string& name) {
    auto& residues = frame.topology().residues();
    std::unordered_set<std::string> string_props;

    for (auto& residue : residues) {
        const auto& prop = residue.get(name);
        if (!prop) {
            continue;
        }
        switch (prop->kind()) {
        case chemfiles::Property::STRING:
            string_props.insert(prop->as_string());
            break;
        default:
            break; // Not a type we support
        }
    }

    return string_props.size();
}

//! Print select residue names and their respective counts
//!
//! \param [in] entry The entry containing the residues of interest
//! \param [in] res_ids The set of residue ids for printing
//! \return a string containing the formatted output.
template <typename Container>
inline std::string print_residue_names(const chemfiles::Frame& entry,
                                       const Container& res_ids) {
    if (res_ids.empty()) {
        return "";
    }

    ResidueNameCount rnc;
    lemon::count::residues(entry, res_ids, rnc);

    std::stringstream ss;
    ss << rnc << "\n";
    return ss.str();
}

} // namespace count

} // namespace lemon

#endif
