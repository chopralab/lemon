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

TEST_CASE("Nitrogen atom typing") {
    chemfiles::Topology topo;
    topo.add_atom(chemfiles::Atom("N"));
    auto type = lemon::xscore::get_xs_type(
        topo, 0, lemon::xscore::create_bond_map(topo.bonds()));
    CHECK(type == lemon::xscore::XS_TYPE_N_D);

    topo.add_atom(chemfiles::Atom("C"));
    topo.add_bond(0, 1, chemfiles::Bond::TRIPLE);
    type = lemon::xscore::get_xs_type(
        topo, 0, lemon::xscore::create_bond_map(topo.bonds()));
    CHECK(type == lemon::xscore::XS_TYPE_N_A);
}

TEST_CASE("Oxygen atom typing") {
    chemfiles::Topology topo;
    topo.add_atom(chemfiles::Atom("O"));
    auto type = lemon::xscore::get_xs_type(
        topo, 0, lemon::xscore::create_bond_map(topo.bonds()));
    CHECK(type == lemon::xscore::XS_TYPE_O_DA);

    topo.add_atom(chemfiles::Atom("C"));
    topo.add_bond(0, 1);
    type = lemon::xscore::get_xs_type(
        topo, 0, lemon::xscore::create_bond_map(topo.bonds()));
    CHECK(type == lemon::xscore::XS_TYPE_O_DA);  // Alcohol

    topo.add_atom(chemfiles::Atom("C"));
    topo.add_bond(0, 2);
    type = lemon::xscore::get_xs_type(
        topo, 0, lemon::xscore::create_bond_map(topo.bonds()));
    CHECK(type == lemon::xscore::XS_TYPE_O_A);  // Ether

    topo.add_atom(chemfiles::Atom("C"));
    topo.add_bond(0, 3);
    type = lemon::xscore::get_xs_type(
        topo, 0, lemon::xscore::create_bond_map(topo.bonds()));
    CHECK(type == lemon::xscore::XS_TYPE_O_P);  // Oddity
}

TEST_CASE("Additional atom typing") {
    chemfiles::Topology topo;
    topo.add_atom(chemfiles::Atom("C"));
    topo.add_atom(chemfiles::Atom("P"));
    topo.add_atom(chemfiles::Atom("F"));
    topo.add_atom(chemfiles::Atom("Br"));
    topo.add_atom(chemfiles::Atom("Cl"));
    topo.add_atom(chemfiles::Atom("I"));
    topo.add_atom(chemfiles::Atom("Mg"));
    topo.add_atom(chemfiles::Atom("Be"));
    topo.add_bond(0, 1);
    topo.add_bond(0, 2);
    topo.add_bond(0, 3);
    topo.add_bond(0, 4);
    topo.add_bond(1, 5);

    auto b_map = lemon::xscore::create_bond_map(topo.bonds());
    CHECK(lemon::xscore::get_xs_type(topo, 1, b_map) ==
          lemon::xscore::XS_TYPE_P_P);
    CHECK(lemon::xscore::get_xs_type(topo, 2, b_map) ==
          lemon::xscore::XS_TYPE_F_H);
    CHECK(lemon::xscore::get_xs_type(topo, 3, b_map) ==
          lemon::xscore::XS_TYPE_Br_H);
    CHECK(lemon::xscore::get_xs_type(topo, 4, b_map) ==
          lemon::xscore::XS_TYPE_Cl_H);
    CHECK(lemon::xscore::get_xs_type(topo, 5, b_map) ==
          lemon::xscore::XS_TYPE_I_H);
    CHECK(lemon::xscore::get_xs_type(topo, 6, b_map) ==
          lemon::xscore::XS_TYPE_Metal_D);
    CHECK(lemon::xscore::get_xs_type(topo, 7, b_map) ==
          lemon::xscore::XS_TYPE_SKIP);
}

TEST_CASE("Scoring") {
    auto traj = chemfiles::Trajectory("files/1AHA.mmtf.gz", 'r');
    auto frame1 = traj.read();

    auto ade = lemon::select::specific_residues(frame1, {"ADE"});
    auto residues = frame1.topology().residues();

    auto small_molecule = residues[*ade.begin()];

    std::set<size_t> proteins;
    for (size_t i = 0; i < residues.size(); ++i) {
        proteins.insert(i);
    }

    lemon::prune::keep_interactions(frame1, ade, proteins, 8.0);
    proteins.erase(*ade.begin());

    auto vscore = lemon::xscore::vina_score(frame1, (*ade.begin()), proteins);
    CHECK(roughly(vscore.g1, 79.23012));
    CHECK(roughly(vscore.g2, 639.67236));
    CHECK(roughly(vscore.hydrogen, 47.031564));
    CHECK(roughly(vscore.hydrophobic, 0.00000));
    CHECK(roughly(vscore.rep, 2.35422));
}

TEST_CASE("Scoring with hydrophobics") {
    auto traj = chemfiles::Trajectory("files/1JD0.mmtf.gz", 'r');
    auto frame1 = traj.read();

    auto azm = lemon::select::specific_residues(frame1, {"AZM"});
    auto residues = frame1.topology().residues();

    auto small_molecule = residues[*azm.begin()];

    std::set<size_t> proteins;
    for (size_t i = 0; i < residues.size(); ++i) {
        proteins.insert(i);
    }

    lemon::prune::keep_interactions(frame1, azm, proteins, 8.0);
    proteins.erase(*azm.begin());

    auto vscore = lemon::xscore::vina_score(frame1, (*azm.begin()), proteins);
    CHECK(roughly(vscore.g1, 52.6882));
    CHECK(roughly(vscore.g2, 766.766, 1e-3));
    CHECK(roughly(vscore.hydrogen, 78.3355));
    CHECK(roughly(vscore.hydrophobic, 17.1033));
    CHECK(roughly(vscore.rep, 4.13918));
}

TEST_CASE("Dry scoring") {
    auto traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame1 = traj.read();

    auto ade = lemon::select::specific_residues(frame1, {"P30"});
    auto residues = frame1.topology().residues();

    auto small_molecule = residues[*ade.begin()];

    std::set<size_t> proteins;
    for (size_t i = 0; i < residues.size(); ++i) {
        proteins.insert(i);
    }

    lemon::prune::identical_residues(frame1, proteins);
    lemon::prune::keep_interactions(frame1, ade, proteins, 8.0);
    proteins.erase(*ade.begin());

    auto vscore = lemon::xscore::vina_score(frame1, (*ade.begin()), proteins);
    CHECK(roughly(vscore.g1, 121.94136));
    CHECK(roughly(vscore.g2, 2103.34799));
    //CHECK(roughly(vscore.hydrogen, 78.3355));
    //CHECK(roughly(vscore.hydrophobic, 17.1033));
    CHECK(roughly(vscore.rep, 14.48645));
}
