#include "benchmarker/select.hpp"
#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Select PSI from 1AAQ") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame = traj.read();

    const auto& res = benchmarker::select_small_molecule(frame);
    CHECK(res.name() == "PSI");
}
