#ifndef LEMON_STRUCTURE_HPP
#define LEMON_STRUCTURE_HPP

#include <algorithm>
#include <set>
#include <limits>

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

inline unsigned long find_operlapping_residues(const chemfiles::Frame& search,
                                        const chemfiles::Frame& native,
                                        std::vector<size_t>& a_search,
                                        std::vector<size_t>& a_native) {
    const auto& search_res = search.topology().residues();
    const auto& native_res = native.topology().residues();

    a_search.reserve(search_res.size());
    a_native.reserve(native_res.size());

    auto n_ali = 0ul;
    for (const auto& s_res : search_res) {
        auto i = find_element_by_name(search, s_res, "CA");
        if (i == search.size()) {
            continue;
        }
        for (const auto& n_res : native_res) {
            if (*s_res.id() == *n_res.id() &&
                s_res.get("chainname") == n_res.get("chainname")) {
                auto j = find_element_by_name(native, n_res, "CA");
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

inline unsigned long count_number_of_atom_names(const chemfiles::Frame& frame,
                                         const std::string& elem) {
    const auto& residues = frame.topology().residues();
    auto c = std::count_if(residues.begin(), residues.end(),
                           [&frame, &elem](const chemfiles::Residue& n_res) {
                               return find_element_by_name(frame, n_res, elem) !=
                                      frame.size();
                           });
    return static_cast<unsigned long>(c);
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

    //! The affine transformation found to maximize the score
    Affine affine;
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
//! \return A strucure with the TMScore, and number of aligned residues
inline TMResult TMscore(const chemfiles::Frame& search,
                        const chemfiles::Frame& native) {
    std::vector<size_t> aligned_search_ids;
    std::vector<size_t> aligned_native_ids;

    // pickup the aligned residues:
    auto n_seq = count_number_of_atom_names(native, "CA");

    // The number of aligned residues
    auto n_ali = find_operlapping_residues(search, native,
                                           aligned_search_ids,
                                           aligned_native_ids);

    if (n_ali == 0) {
        return {0.0, 0.0, 0, Affine{{0.0, 0.0, 0.0}, Matrix3D::unit()}};
    }

    // parameters:
    auto d0 = 1.24 * std::pow(n_seq - 15, 1.0 / 3.0) - 1.8;
    d0 = std::max(0.5, d0);

    auto d0_search = clamp(d0, 4.5, 8.0);

    // iterative parameters
    auto max_iter = 20ul;      // maximum number of iterations
    auto local_init_max = 6ul; // maximum number of local starting points

    auto local_init = std::vector<size_t>();
    local_init.reserve(local_init_max);

    auto local_init_min = std::min(n_ali, 4ul);

    // Setup up the number of local structure iterations in local_init
    // The last value of this array must be less than or equal to local_init_min
    auto n_init = 0ul; // Total number of local iterations. [0, local_init_max]
    while (n_init < local_init_max - 1) {
        auto temp = std::pow(2.0, -1.0 * static_cast<double>(n_init));
        temp *= static_cast<double>(n_ali);
        local_init.push_back(static_cast<size_t>(temp));
        if (local_init[n_init] <= local_init_min) {
            local_init[n_init] = local_init_min;
            break;
        }

        ++n_init;
    }

    // Check to see if the loop exited early
    if (n_init == local_init_min - 1) {
        local_init.push_back(local_init_min);
    }

    // find the maximum score starting from local structures superposition
    const auto& search_pos = search.positions();
    const auto& native_pos = native.positions();

    // Calculate the TMscore of an `align_pos` to the `native_pos`.
    auto score_fun = [&](double d, const Coordinates& align_pos) {
        auto score_sum = 0.0;
        auto cut_scored_ids = std::vector<size_t>(n_ali);

        // The number of well aligned residues. Updated by score_fun
        auto n_cut = 0ul;

        do {
            n_cut = 0; // number of residue-pairs dis<d, for iteration
            score_sum = 0.0;
            for (size_t k = 0; k < n_ali; ++k) {
                auto j = aligned_native_ids[k];

                auto dis = (align_pos[k] - native_pos[j]).norm();
                if (dis < d) {
                    // [0,n_ali - 1], mark the residue-pairs in dis<d
                    cut_scored_ids[n_cut++] = k;
                }
                score_sum += 1 / (1 + std::pow(dis / d0, 2));
            }
            d += 0.5;
        } while (n_cut < 3 && n_ali > 3);

        cut_scored_ids.resize(n_cut);

        return std::make_pair(score_sum / static_cast<double>(n_seq), cut_scored_ids);
    };

    // Aligns a given number of positions to native positions and retuns
    // the result.
    // If `start` is zero, the `cut_scored_ids` vector is used to find aligned
    // residues, otherwise it is calculated from `start` and the current loop index
    auto rotate_substructure = [&](size_t start, size_t length,
                                   const std::vector<size_t>& cut_ids) {

        // Which residues have been aligned
        auto aligned_ids = std::vector<size_t>(length);

        auto rotated_search_pos = Coordinates(n_ali);

        for (auto k = 0ul; k < n_ali; ++k) {
            auto i = aligned_search_ids[k];
            rotated_search_pos[k] = search_pos[i];
        }

        auto search_align_pos = Coordinates(length);
        auto native_align_pos = Coordinates(length);

        for (size_t i = 0; i < length; ++i) {
             // The current residue to copy into the alignment arrays
            auto k = start == 0 ? cut_ids[i] : start - 1 + i;

            search_align_pos[i] = search_pos[aligned_search_ids[k]];
            native_align_pos[i] = native_pos[aligned_native_ids[k]];

            aligned_ids[i] = k;
        }

        auto affine = kabsch(std::move(native_align_pos),
                             std::move(search_align_pos));

        lemon::align(rotated_search_pos, affine);

        return std::make_tuple(rotated_search_pos, aligned_ids, affine);
    };

    auto score_max = std::numeric_limits<double>::min();

    auto best_affine = Affine{{0.0, 0.0, 0.0}, Matrix3D::unit()};

    for (auto local_start : local_init) {

        for (auto iL = 1ul; iL <= n_ali - local_start + 1; ++iL) {

            double score;
            auto cut_ids = std::vector<size_t>();
            auto aligned = std::vector<size_t>();
            auto rotated = Coordinates();
            auto affine = Affine{{0.0, 0.0, 0.0}, Matrix3D::unit()};

            // Rotate based on current alignment
            std::tie(rotated, aligned, affine) = rotate_substructure(iL, local_start, cut_ids);

            // Get the current score and update `cut_ids`
            std::tie(score, cut_ids) = score_fun(d0_search - 1, rotated);

            if (score > score_max) {
                score_max = score;
                best_affine = std::move(affine);
            }

            // iterations for extending the local search
            for (auto it = 0ul; it <= max_iter; ++it) {
                std::tie(rotated, aligned, affine) = rotate_substructure(0, cut_ids.size(), cut_ids);

                // get scores, n_cut+cut_scored_ids(i) for iteration
                std::tie(score, cut_ids) = score_fun(d0_search + 1, rotated);

                if (score > score_max) {
                    score_max = score;
                    best_affine = std::move(affine);
                }

                if (cut_ids == aligned) {
                    break;
                }
            }
        }
    }

    return {score_max, 0.0, n_ali, best_affine};
}

} // namespace tmalign

} // namespace lemon

#endif
