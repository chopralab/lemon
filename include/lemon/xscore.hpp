#ifndef LEMON_SCORE_HPP
#define LEMON_SCORE_HPP

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#include <chemfiles/Frame.hpp>
LEMON_EXTERNAL_FILE_POP

#include <cmath>
#include <set>
#include <unordered_map>

namespace lemon {

namespace xscore {

enum class XS_TYPE {
    C_H = 0,
    C_P = 1,
    N_P = 2,
    N_D = 3,
    N_A = 4,
    N_DA = 5,
    O_P = 6, // Unused???
    O_D = 7,
    O_A = 8,
    O_DA = 9,
    S_P = 10,
    P_P = 11,
    F_H = 12,
    Cl_H = 13,
    Br_H = 14,
    I_H = 15,
    Metal_D = 16,
    SKIP = 17, // Skip for calculations
};

const double vdw_radii[] = {
    1.9, // C_H
    1.9, // C_P
    1.8, // N_P
    1.8, // N_D
    1.8, // N_A
    1.8, // N_DA
    1.7, // O_P
    1.7, // O_D
    1.7, // O_A
    1.7, // O_DA
    2.0, // S_P
    2.1, // P_P
    1.5, // F_H
    1.8, // Cl_H
    2.0, // Br_H
    2.2, // I_H
    1.2, // Metal_D
    0.0, // Skip
};

using lemon::xscore::XS_TYPE;

inline bool is_hydrophobic(XS_TYPE xs) {
    switch (xs) {
    case XS_TYPE::C_H:
    case XS_TYPE::F_H:
    case XS_TYPE::Cl_H:
    case XS_TYPE::Br_H:
    case XS_TYPE::I_H:
        return true;
    default:
        return false;
    }
}

inline bool is_acceptor(XS_TYPE xs) {
    switch (xs) {
    case XS_TYPE::N_A:
    case XS_TYPE::N_DA:
    case XS_TYPE::O_A:
    case XS_TYPE::O_DA:
        return true;
    default:
        return false;
    }
}

inline bool is_donor(XS_TYPE xs) {
    switch (xs) {
    case XS_TYPE::N_D:
    case XS_TYPE::N_DA:
    case XS_TYPE::O_D:
    case XS_TYPE::O_DA:
    case XS_TYPE::Metal_D:
        return true;
    default:
        return false;
    }
}

inline bool donor_acceptor(XS_TYPE t1, XS_TYPE t2) {
    return is_donor(t1) && is_acceptor(t2);
}

inline bool h_bond_possible(XS_TYPE t1, XS_TYPE t2) {
    return donor_acceptor(t1, t2) || donor_acceptor(t2, t1);
}

inline double optimal_distance(XS_TYPE xs_t1, XS_TYPE xs_t2) {
    return vdw_radii[static_cast<size_t>(xs_t1)] +
           vdw_radii[static_cast<size_t>(xs_t2)];
}

//! An object which represents the components of XScore/Vina's scoring function
//!
//! The XScore/Vina scoring function is divided into five components:
//! two gaussian functions, hydrophobic interactions, repulsive forces, and
//! hydrogen bonding.
struct VinaScore {
    double g1 = 0;          //!< The small gaussian term
    double g2 = 0;          //!< The large gaussian term
    double rep = 0;         //!< The repulsive term
    double hydrophobic = 0; //!< The hydrophobic term
    double hydrogen = 0;    //!< The hydrogen bonding term
};

using BondMap = std::unordered_multimap<size_t, size_t>;

// Easy case, just check if carbon is bound to a heteroatom
inline XS_TYPE
get_c_xs_type(const chemfiles::Topology& topo, size_t j,
              const BondMap& bond_map) {
    auto range = bond_map.equal_range(j);
    auto& bonds = topo.bonds();

    for (auto k = range.first; k != range.second; ++k) {
        auto bond = bonds[k->second];
        auto neighbor = (j == bond[0]) ? bond[1] : bond[0];
        auto type = *(topo[neighbor].atomic_number());
        if (!(type == 6 || type == 1)) {
            return XS_TYPE::C_P;
        }
    }

    return XS_TYPE::C_H;
}

inline bool is_in_strong_resonance(
    const chemfiles::Topology& topo, size_t j,
    const BondMap& bond_map) {
    auto range = bond_map.equal_range(j);
    auto& bond_order = topo.bond_orders();
    auto& bonds = topo.bonds();

    for (auto k = range.first; k != range.second; ++k) {
        auto bond = bonds[k->second];
        auto neighbor = (j == bond[0]) ? bond[1] : bond[0];
        auto type = *(topo[neighbor].atomic_number());
        if ((type == 8 || type == 7) && bond_order[k->second] == 2) {
            return true;
        }
    }

    return false;
}

// Hard case, there's three bond counts to get, which in turn involve checks
inline XS_TYPE
get_n_xs_type(const chemfiles::Topology& topo, size_t j,
              const BondMap& bond_map) {
    auto range = bond_map.equal_range(j);
    auto& bond_order = topo.bond_orders();
    auto& bonds = topo.bonds();

    if (bond_map.count(j) == 0) { // Typically Ammonia, assume protonated
        return XS_TYPE::N_D;
    } else if (bond_map.count(j) == 1) { // Nitrile, monoamine, etc
        auto bond = bond_map.find(j);
        if (bond_order[bond->second] == 3) {
            return XS_TYPE::N_A; // It can only accept - Nitrile
        }
        
        if (bond_order[bond->second] == 2) {
            return XS_TYPE::N_D; // Imine, or guanidine
        }

        return XS_TYPE::N_D; // AutoDOCK labels amines as Donor only
    } else { // A linker, amide, guanidine, amine of some sort, etc,
        size_t sum_of_orders = 0;
        size_t sum_of_hydrogens = 0;
        size_t is_withdraw = 0;
        for (auto k = range.first; k != range.second; ++k) {
            auto bond = bonds[k->second];
            auto neighbor = (j == bond[0]) ? bond[1] : bond[0];
            auto type = *(topo[neighbor].atomic_number());
            sum_of_hydrogens += static_cast<size_t>(type == 1); // Sum explicit hydrogens
            sum_of_orders += static_cast<size_t>(bond_order[k->second]);
            is_withdraw +=
                static_cast<size_t>(is_in_strong_resonance(topo, neighbor, bond_map));
        }

        if (sum_of_hydrogens >= 2) {
            return XS_TYPE::N_D;
        }

        if (sum_of_orders == 2 && bond_map.count(j) == 2) {
            // Linking Amine/Amide or heterocycle
            return is_withdraw != 0U ? XS_TYPE::N_P : XS_TYPE::N_D;
        }

        if (sum_of_orders == 3 && bond_map.count(j) == 2) {
            return XS_TYPE::N_A;
        }

        if (sum_of_orders == 3 && is_withdraw == 0U) {
            return XS_TYPE::N_A; // Linking SP2 nitrogen, aromatic, etc
        }
        // Fall through to Polar nitrogen, ie C=N=C, next to withdrawing group
    }
    return XS_TYPE::N_P;
}

// Medium case, two bond cases - but no additional checks.
inline XS_TYPE
get_o_xs_type(const chemfiles::Topology& topo, size_t j,
              const BondMap& bond_map) {
    auto range = bond_map.equal_range(j);
    auto& bond_order = topo.bond_orders();
    auto& bonds = topo.bonds();

    if (bond_map.count(j) == 0) { // Typically water
        return XS_TYPE::O_DA;
    } else if (bond_map.count(j) == 1) { // Carbonyl or alcohol
        auto bond = bond_map.find(j);
        auto bond2 = bonds[bond->second];
        auto neighbor = (j == bond2[0]) ? bond2[1] : bond2[0];

        if (bond_order[bond->second] == 2) {
            return XS_TYPE::O_A; // Carbonyl
        }

        if (is_in_strong_resonance(topo, neighbor, bond_map)) {
            return XS_TYPE::O_A; // Carboxylic acid, phosphate, sulfate
        }

        return XS_TYPE::O_DA; // Alcohol
    } else if (bond_map.count(j) == 2) { // alcohol or ether
        for (auto k = range.first; k != range.second; ++k) {
            auto bond = bonds[k->second];
            auto neighbor = (j == bond[0]) ? bond[1] : bond[0];
            auto type = *(topo[neighbor].atomic_number());
            if (type == 1) {
                return XS_TYPE::O_DA; // Alcohol with explict H
            }
        }
        return XS_TYPE::O_A; // ether, phosphodiester, etc
    }
    // Oddity, not possible in RCSB
    return XS_TYPE::O_P;
}

inline XS_TYPE
get_xs_type(const chemfiles::Topology& topo, size_t j,
            const BondMap& bond_map) {
    auto atomic_id = *(topo[j].atomic_number());

    switch (atomic_id) {
    // Handle 'easy' cases
    case 9:  // NOLINT Fluorine
        return XS_TYPE::F_H;
    case 15: // NOLINT Phosphorous
        return XS_TYPE::P_P;
    case 16: // NOLINT Sulfur
        return XS_TYPE::S_P;
    case 17: // NOLINT Chlorine
        return XS_TYPE::Cl_H;
    case 35: // NOLINT Bromine
        return XS_TYPE::Br_H;
    case 53: // NOLINT Iodine
        return XS_TYPE::I_H;
    case 12: // NOLINT Mg
    case 20: // NOLINT Ca
    case 25: // NOLINT Mn
    case 26: // NOLINT Fe
    case 30: // NOLINT Zn
        return XS_TYPE::Metal_D;
    case 6:
        return get_c_xs_type(topo, j, bond_map);
    case 7:
        return get_n_xs_type(topo, j, bond_map);
    case 8:
        return get_o_xs_type(topo, j, bond_map);
    }
    return XS_TYPE::SKIP;
}

inline double gaussian(double offset, double width, double r) {
    double exponent = (r - offset) / width;
    return std::exp(-exponent * exponent);
}

inline double repulsion(double offset, double r) {
    r -= offset;
    if (r > 0.0) {
        return 0.0;
    }

    return r * r;
}

inline double slope_step(double end, double intercept, double r) {
    if(intercept < end) {
        if(r <= intercept) {
            return 0;
        }
        if(r >= end) {
            return 1;
        }
    } else {
        if(r >= intercept) {
            return 0;
        }
        if(r <= end) {
            return 1;
        }
    }
    return (r - intercept) / (end - intercept);
}

// Move to a topology.hpp file?
inline BondMap create_bond_map(const std::vector<chemfiles::Bond>& bonds) {
    BondMap bond_map;

    // Map bonds
    for (size_t i = 0; i < bonds.size(); ++i) {
        auto& bond = bonds[i];
        bond_map.insert({bond[0], i});
        bond_map.insert({bond[1], i});
    }

    return bond_map;
}

//! Default maximum for interactions
constexpr auto DEFAULT_INTERACTION_DISTANCE = 8.0;

//! XScore is a 'docking' scoring function used to evaluate compound-protein
//! interactions
//!
//! The docking program AutoDOCK Vina utilizes a modified version of the XScore
//! scoring function to evaluate the fit of a compound-protein interaction.
//! Since the original XScore program is not open source, we've included a copy
//! of the Vina version of this scoring function.
//! \param [in] frame The frame for which the ligand-protein score will be
//!  calculated
//! \param [in] ligid The residue ID for the ligand in the entry
//! \param [in] recid The residue IDs for the protein in the entry
//! \param [in] cutoff The interaction distance cutoff between ligand and
//!  protein
//! \return The five components of Vina/XScore's scoring function.
template <typename Container>
inline VinaScore vina_score(const chemfiles::Frame& frame, size_t ligid,
                            Container recid,
                            double cutoff = DEFAULT_INTERACTION_DISTANCE) {
    const auto& topo = frame.topology();
    const auto& residues = topo.residues();
    auto bond_map = create_bond_map(topo.bonds());

    // Not memory efficient, but simple to implement
    std::vector<XS_TYPE> xs_types(frame.size());
    for (auto i : recid) {
        for (auto j : residues[i]) {
            xs_types[j] = get_xs_type(topo, j, bond_map);
        }
    }

    VinaScore X_Score;
    auto& small_molecule = residues[ligid];
    for (auto i : small_molecule) {
        xs_types[i] = get_xs_type(topo, i, bond_map);
        if (xs_types[i] == XS_TYPE::SKIP) {
            continue;
        }

        for (auto k : recid) {
            for (auto j : residues[k]) {
                if (xs_types[j] == XS_TYPE::SKIP) {
                    continue;
                }

                auto r = frame.distance(i, j);
                if (r > cutoff) {
                    continue;
                }

                r -= optimal_distance(xs_types[i], xs_types[j]);
                X_Score.g1 += gaussian(0, 0.5, r);
                X_Score.g2 += gaussian(3, 2.0, r);

                X_Score.rep += repulsion(0, r);
                if (is_hydrophobic(xs_types[i]) &&
                    is_hydrophobic(xs_types[j])) {
                    X_Score.hydrophobic += slope_step(0.5, 1.5, r); // NOLINT
                }
                if (h_bond_possible(xs_types[i], xs_types[j])) {
                    X_Score.hydrogen += slope_step(-0.7, 0, r); // NOLINT
                }
            }
        }
    }

    return X_Score;
}

} // namespace xscore

} // namespace lemon

#endif
