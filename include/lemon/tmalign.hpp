#ifndef LEMON_STRUCTURE_HPP
#define LEMON_STRUCTURE_HPP

#include <algorithm>
#include <set>

#include "chemfiles/Frame.hpp"
#include "lemon/matrix.hpp"

namespace lemon {

namespace tmalign {

inline size_t find_element_by_name(const chemfiles::Frame& frame,
                                   const chemfiles::Residue& res,
                                   const std::string& elem) {
    for (auto i : res) {
        if (frame[i].name() == elem) {
            return i;
        }
    }

    return frame.size();
}

inline size_t find_operlapping_residues(const chemfiles::Frame& search,
                                        const chemfiles::Frame& native,
                                        std::vector<size_t>& a_search,
                                        std::vector<size_t>& a_native) {
    const auto& search_res = search.topology().residues();
    const auto& native_res = native.topology().residues();

    a_search.reserve(search_res.size());
    a_native.reserve(native_res.size());

    size_t n_ali = 0;
    for (const auto& sres : search_res) {
        auto i = find_element_by_name(search, sres, "CA");
        if (i == search.size()) {
            continue;
        }
        for (const auto& nres : native_res) {
            if (*sres.id() == *nres.id() &&
                sres.get("chainname") == nres.get("chainname")) {
                auto j = find_element_by_name(native, nres, "CA");
                if (j == native.size()) {
                    continue;
                }

                a_search.emplace_back(i);
                a_native.emplace_back(j);
                n_ali++;
                break;
            }
        }
    }

    return n_ali;
}

inline size_t count_number_of_atom_names(const chemfiles::Frame& frame,
                                         const std::string& elem) {
    const auto& residues = frame.topology().residues();
    auto c = std::count_if(residues.begin(), residues.end(),
                           [&frame, &elem](const chemfiles::Residue& nres) {
                               return find_element_by_name(frame, nres, elem) !=
                                      frame.size();
                           });
    return static_cast<size_t>(c);
}

//! Structure representing the TMScore
//!
//!
struct TMResult {
    //! The TMscore
    double score;

    //! The RMSD of the aligned proteins
    double rmsd;

    //! Number of aligned residues
    size_t aligned;
};

//! TMalign is an algorithm used to align protein chains in 3D space.
//!
//! Ported from
//! https://zhanglab.ccmb.med.umich.edu/TM-score/TMscore_subroutine.f Original
//! reference: Yang Zhang, Jeffrey Skolnick, Proteins 2004 57:702-10.
//!
//! Original License:
//! Permission to use, copy, modify, and distribute this program for
//! any purpose, with or without fee, is hereby granted, provided that
//! the notices on the head, the reference information, and this
//! copyright notice appear in all copies or substantial portions of
//! the Software. It is provided "as is" without express or implied
//! warranty.
//! \param [in] search The frame that is being aligned to native.
//! \param [in] native The 'native' chain that the search chain is aligned to.
//! \param [out] rot The aligned version of the search frame.
//! \param [in] align Should the search frame to aligned to native.
//! \return A strucure with theTMScore, RMSD after alignment, and number of
//! aligned residues
inline TMResult TMscore(const chemfiles::Frame& search,
                        const chemfiles::Frame& native) {
    std::vector<size_t> a_search;
    std::vector<size_t> a_native;

    // pickup the aligned residues:
    auto nseqB = count_number_of_atom_names(native, "CA");
    auto n_ali = find_operlapping_residues(search, native, a_search, a_native);

    if (n_ali == 0) {
        return {0.0, 0.0, 0};
    }

    // parameters:
    // d0------------->
    double d0 = 1.24 * std::pow(nseqB - 15, 1.0 / 3.0) - 1.8;
    if (d0 < 0.5)
        d0 = 0.5;

    // d0_search ----->
    double d0_search = d0;
    if (d0_search > 8)
        d0_search = 8;
    if (d0_search < 4.5)
        d0_search = 4.5;

    // iterative parameters ----->
    size_t n_it = 20;      // maximum number of iterations
    size_t n_init_max = 6; // maximum number of L_init
    size_t L_ini_min = 4;
    std::vector<size_t> L_ini(n_init_max);
    if (n_ali < 4)
        L_ini_min = n_ali;

    size_t n_init;
    for (n_init = 0; n_init < n_init_max; ++n_init) {
        if (n_init_max == n_init_max - 1) {
            L_ini[n_init] = L_ini_min;
            break;
        }

        auto temp = std::pow(2.0, -1.0 * static_cast<double>(n_init));
        temp *= static_cast<double>(n_ali);
        L_ini[n_init] = static_cast<size_t>(temp);
        if (L_ini[n_init] <= L_ini_min) {
            L_ini[n_init] = L_ini_min;
            break;
        }
    }

    // find the maximum score starting from local structures superposition
    const auto& search_pos = search.positions();
    const auto& native_pos = native.positions();
    std::vector<double> w(n_ali, 1.0);

    size_t n_cut;
    std::vector<size_t> i_ali(n_ali);
    std::vector<size_t> k_ali(n_ali);
    std::vector<size_t> k_ali0(n_ali);

    Coordinates rot;
    rot.resize(search.size());

    for (size_t k = 0; k < n_ali; ++k) {
        auto i = a_search[k];
        rot[i] = search_pos[i];
    }

    auto score_fun = [&](double d) -> double {
        double score_sum;
        do {
            n_cut = 0; // number of residue-pairs dis<d, for iteration
            score_sum = 0.0;
            for (size_t k = 0; k < n_ali; ++k) {
                auto i = a_search[k]; // [1,nseqA] reoder number of structureA
                auto j = a_native[k]; // [1,nseqB]

                auto diff = rot[i] - native_pos[j];
                double dis = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] +
                                       diff[2] * diff[2]);
                if (dis < d) {
                    // [0,n_ali - 1], mark the residue-pairs in dis<d
                    i_ali[n_cut++] = k;
                }
                score_sum += 1 / (1 + std::pow(dis / d0, 2));
            }
            d += 0.5;
        } while (n_cut < 3 && n_ali > 3);
        return score_sum / static_cast<double>(nseqB);
    };

    double score_max = -1;
    size_t ka0 = 0;

    for (size_t i_init = 0; i_init < n_init; ++i_init) {
        auto L_init = L_ini[i_init];
        auto iL_max = n_ali - L_init + 1;

        std::vector<chemfiles::Vector3D> r_1(n_ali);
        std::vector<chemfiles::Vector3D> r_2(n_ali);

        for (size_t iL = 0; iL < iL_max; ++iL) { // on residues, [0, nseqA - 1]
            size_t ka = 0;

            for (size_t i = 0; i < L_init; ++i) {
                auto k = iL + i; // [0,n_ali - 1] common aligned

                r_1[i] = search_pos[a_search[k]];
                r_2[i] = native_pos[a_native[k]];

                k_ali[ka++] = k;
            }

            r_1.resize(L_init);
            r_2.resize(L_init);

            auto affine = kabsch(r_2, r_1);

            for (size_t k = 0; k < n_ali; ++k) {
                auto j = a_search[k];
                rot[j] = search_pos[j];
            }

            lemon::align(rot, affine);

            // init, get scores, n_cut+i_ali(i) for iteration
            auto score = score_fun(d0_search - 1);
            if (score_max < score) {
                score_max = score;
                ka0 = ka;
                for (size_t i = 0; i < ka0; ++i) {
                    k_ali0[i] = k_ali[i];
                }
            }

            //   iteration for extending ---------------------------------->
            for (size_t it = 1; it <= n_it; ++it) {
                ka = 0;
                r_1.resize(n_cut);
                r_2.resize(n_cut);
                for (size_t i = 0; i < n_cut; ++i) {
                    auto m = i_ali[i]; // [0,n_ali - 1]
                    r_1[i] = search_pos[a_search[m]];
                    r_2[i] = native_pos[a_native[m]];
                    k_ali[ka++] = m;
                }

                // u rotate r_1 to r_2
                affine = kabsch(r_2, r_1);
                for (size_t k = 0; k < n_ali; ++k) {
                    auto j = a_search[k];
                    rot[j] = search_pos[j];
                }

                lemon::align(rot, affine);

                // get scores, n_cut+i_ali(i) for iteration
                score = score_fun(d0_search + 1);

                if (score_max < score) {
                    score_max = score;
                    ka0 = ka;
                    for (size_t i = 0; i <= ka; ++i) {
                        k_ali0[i] = k_ali[i];
                    }
                }
                if (it == n_it)
                    break;
                if (n_cut == ka) {
                    size_t neq = 0;
                    for (size_t i = 0; i < n_cut; ++i) {
                        if (i_ali[i] == k_ali[i])
                            ++neq;
                    }
                    if (n_cut == neq)
                        break;
                }
            }
        } //! for shift
    }     //! for initial length, L_ali/M

    return {score_max, 0.0, n_ali};
}

} // namespace tmalign

} // namespace lemon

#endif
