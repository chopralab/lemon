#include "lemon/score.hpp"
#include <chemfiles.hpp>
#include "lemon/prune.hpp"
#include "lemon/select.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

static bool roughly(double a, double b, double tol = 1e-4) {
    return std::fabs(a - b) < tol;
}

TEST_CASE("Gaussian Terms") {
    double r =
        6.96584 - lemon::xscore::optimal_distance(lemon::xscore::XS_TYPE_N_A,
                                                  lemon::xscore::XS_TYPE_O_A);
    CHECK(roughly(lemon::xscore::gaussian(3, 2, r), 0.947194));
}

TEST_CASE("Scoring") {
    auto traj = chemfiles::Trajectory("files/1AHA.mmtf.gz", 'r');
    auto frame1 = traj.read();

    std::unordered_multimap<size_t, size_t> bond_map;

    // Map bonds
    auto& bonds = frame1.topology().bonds();
    for (size_t i = 0; i < bonds.size(); ++i) {
        auto& bond = bonds[i];
        bond_map.insert({bond[0], i});
        bond_map.insert({bond[1], i});
    }

    auto ade = lemon::select_specific_residues(frame1, {"ADE"});
    auto residues = frame1.topology().residues();

    auto small_molecule = residues[*ade.begin()];
    std::vector<size_t> xs_types(small_molecule.size());
    size_t count = 0;
    for (auto i : small_molecule) {
        xs_types[count] =
            lemon::xscore::get_xs_type(frame1, i, bond_map);
        std::cout << i << " " << xs_types[count] << std::endl;
        ++count;
    }

    std::set<size_t> proteins;
    for (size_t i = 0; i < residues.size(); ++i) {
        proteins.insert(i);
    }

    lemon::keep_interactions(frame1, ade, proteins, 8.0);
    proteins.erase(*ade.begin());

    auto vscore = lemon::xscore::vina_score(frame1, (*ade.begin()), proteins);
    std::cout << vscore.g1 << " " << vscore.g2 << " " << vscore.rep << " "
              << vscore.hydrophobic << " " << vscore.hydrogen << std::endl;
}