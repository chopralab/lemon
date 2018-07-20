#include "lemon/prune.hpp"
#include <chemfiles.hpp>
#include "lemon/select.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Select and prune P30 from 4XUF") {
    auto traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame = traj.read();

    auto res = lemon::select_small_molecule(frame);
    CHECK(res.size() == 2);

    lemon::remove_identical_residues(frame, res);
    CHECK(res.size() == 1);
}

TEST_CASE("Select and prune 1PG/HEM from 1D7D") {
    auto traj = chemfiles::Trajectory("files/1D7D.mmtf.gz", 'r');
    auto frame = traj.read();

    auto res = lemon::select_small_molecule(frame);
    CHECK(res.size() == 3);  // Two hemes and 1PG

    lemon::remove_cofactors(frame, res, lemon::linear_molecules);
    CHECK(res.size() == 2);

    lemon::remove_cofactors(frame, res, lemon::common_cofactors);
    CHECK(res.size() == 0);
}

TEST_CASE("Remove non-nucleic acid interactions") {
    auto traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame = traj.read();

    auto res = lemon::select_small_molecule(frame);
    auto nas = lemon::select_nucleic_acids(frame);
    lemon::find_interactions(frame, res, nas);
    CHECK(res.size() == 0);  // It got removed

    traj = chemfiles::Trajectory("files/entry_10/1/0/100D.mmtf.gz", 'r');
    frame = traj.read();

    res = lemon::select_small_molecule(frame);
    nas = lemon::select_nucleic_acids(frame);
    lemon::find_interactions(frame, res, nas);
    CHECK(res.size() == 1);  // Not removed
}

TEST_CASE("Remove non-metal interactions") {
    auto traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame = traj.read();

    auto res = lemon::select_small_molecule(frame);
    auto metals = lemon::select_metal_ions(frame);
    lemon::find_interactions(frame, res, metals);
    CHECK(res.size() == 0);  // It got removed

    traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    frame = traj.read();

    res = lemon::select_small_molecule(frame);
    metals = lemon::select_metal_ions(frame);
    lemon::find_interactions(frame, res, metals);
    CHECK(res.size() == 1);  // Not removed
}
