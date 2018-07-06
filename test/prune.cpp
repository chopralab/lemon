#include "benchmarker/prune.hpp"
#include "benchmarker/select.hpp"
#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Select and prune P30 from 4XUF") {
    auto traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame = traj.read();

    auto res = benchmarker::select_small_molecule(frame);
    CHECK(res.size() == 2);

    benchmarker::remove_identical_residues(frame, res);
    CHECK(res.size() == 1);
}
