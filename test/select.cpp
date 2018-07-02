#include "benchmarker/select.hpp"
#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Select PSI from 1AAQ") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame = traj.read();

    const auto& res = benchmarker::select_small_molecule(frame);
    CHECK(res.size() == 1);

    size_t id = *res.begin();
    CHECK(frame.topology().residues()[id].name() == "PSI");
}

TEST_CASE("Select FE from 2WTL") {
    auto traj = chemfiles::Trajectory("files/2WTL.mmtf.gz", 'r');
    auto frame = traj.read();

    const auto& res = benchmarker::select_metal_ions(frame);
    CHECK(res.size() == 12);
}
