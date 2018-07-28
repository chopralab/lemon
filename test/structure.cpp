#include "lemon/structure.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/Trajectory.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("TMScore") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    double TMscore, rmsd;
    size_t aligned_res; 
    std::tie(TMscore, rmsd, aligned_res) = lemon::TMscore(frame1, frame2, rot);
}