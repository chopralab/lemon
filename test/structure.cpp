#include "lemon/structure.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/Trajectory.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

static bool roughly(double a, double b, double tol = 1e-6) {
    return std::fabs(a - b) < tol;
}

TEST_CASE("TMScore for 1OQ5") {
    auto traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    double TMscore, rmsd;
    size_t aligned_res; 
    std::tie(TMscore, rmsd, aligned_res) = lemon::TMscore(frame1, frame2, rot);
    CHECK(roughly(TMscore, 1));
    CHECK(roughly(rmsd, 0));
    CHECK(aligned_res == 256);
}

TEST_CASE("TMScore for 1AAQ") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    double TMscore, rmsd;
    size_t aligned_res; 
    std::tie(TMscore, rmsd, aligned_res) = lemon::TMscore(frame1, frame2, rot);
    CHECK(roughly(TMscore, 1));
    CHECK(roughly(rmsd,0));
    CHECK(aligned_res == 199);
}

TEST_CASE("TMScore for 1OQ5/1AAQ") {
    auto traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    double TMscore, rmsd;
    size_t aligned_res; 
    std::tie(TMscore, rmsd, aligned_res) = lemon::TMscore(frame1, frame2, rot);
    CHECK(roughly(TMscore, 0.100864));
    CHECK(roughly(rmsd, 16.0814691307));
    CHECK(aligned_res == 96);
}

TEST_CASE("TMScore for 1AAQ/1OQ5") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    double TMscore, rmsd;
    size_t aligned_res; 
    std::tie(TMscore, rmsd, aligned_res) = lemon::TMscore(frame1, frame2, rot, true);
    CHECK(roughly(TMscore, 0.0863468016));
    CHECK(roughly(rmsd, 16.0814691307));
    CHECK(aligned_res == 96);
}


TEST_CASE("TMScore for 1AAQ/4XUF") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    double TMscore, rmsd;
    size_t aligned_res; 
    std::tie(TMscore, rmsd, aligned_res) = lemon::TMscore(frame1, frame2, rot);
    CHECK(roughly(TMscore,0.0));
    CHECK(roughly(rmsd, 0));
    CHECK(aligned_res == 0);
}
