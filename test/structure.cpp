#include "lemon/structure.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/Trajectory.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

static bool roughly(double a, double b, double tol = 1e-4) {
    return std::fabs(a - b) < tol;
}

TEST_CASE("tmalign::TMscore for 1OQ5") {
    auto traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    double tmscore, rmsd;
    size_t aligned_res;
    auto tm =
        lemon::tmalign::TMscore(frame1, frame2, rot);
    CHECK(roughly(tm.score, 1));
    CHECK(roughly(tm.rmsd, 0));
    CHECK(tm.aligned == 256);
}

TEST_CASE("tmalign::TMscore for 1AAQ") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    auto tm =
        lemon::tmalign::TMscore(frame1, frame2, rot);
    CHECK(roughly(tm.score, 1));
    CHECK(roughly(tm.rmsd, 0));
    CHECK(tm.aligned == 199);
}

TEST_CASE("tmalign::TMscore for 1OQ5/1AAQ") {
    auto traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    auto tm =
        lemon::tmalign::TMscore(frame1, frame2, rot);
    CHECK(roughly(tm.score, 0.100864));
    CHECK(roughly(tm.rmsd, 16.0814691307));
    CHECK(tm.aligned == 96);
}

TEST_CASE("tmalign::TMscore for 1AAQ/1OQ5") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    auto tm =
        lemon::tmalign::TMscore(frame1, frame2, rot, true);
    CHECK(roughly(tm.score, 0.0863468016));
    CHECK(roughly(tm.rmsd, 16.0814691307));
    CHECK(tm.aligned == 96);
}

TEST_CASE("tmalign::TMscore for 1AAQ/4XUF") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    auto tm =
        lemon::tmalign::TMscore(frame1, frame2, rot);
    CHECK(roughly(tm.score, 0.0));
    CHECK(roughly(tm.rmsd, 0));
    CHECK(tm.aligned == 0);
}

TEST_CASE("tmalign::TMscore for 1AAQ/1YT9") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1YT9.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<chemfiles::Vector3D> rot;
    auto tm =
        lemon::tmalign::TMscore(frame1, frame2, rot);
    CHECK(roughly(tm.score, 0.9762572135));
    CHECK(roughly(tm.rmsd, 0.8351581378));
    CHECK(tm.aligned == 198);
}

TEST_CASE("Kabsch") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1YT9.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<size_t> a_search;
    std::vector<size_t> a_native;

    lemon::tmalign::find_operlapping_residues(frame1, frame2, a_search, a_native);

    const auto& a = frame1.positions();
    const auto& b = frame2.positions();

    std::vector<chemfiles::Vector3D> r_1(a_search.size());
    std::vector<chemfiles::Vector3D> r_2(a_search.size());

    for (size_t i = 0; i < a_search.size(); ++i) {
        r_1[i] = a[a_search[i]];
        r_2[i] = b[a_native[i]];
    }

    int ier;
    std::vector<double> w(a_search.size(), 1.0);
    chemfiles::Vector3D t;
    chemfiles::Matrix3D u = chemfiles::Matrix3D::zero();
    lemon::tmalign::kabsch(w, r_1, r_2, a_search.size(), u, t, ier, 2.0);
    lemon::tmalign::kabsch(w, r_1, r_2, a_search.size(), u, t, ier, 1.0);
}

TEST_CASE("Kabsch 2") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1YT9.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<size_t> a_search;
    std::vector<size_t> a_native;

    lemon::tmalign::find_operlapping_residues(frame1, frame2, a_search, a_native);

    const auto& a = frame1.positions();
    const auto& b = frame2.positions();

    std::vector<chemfiles::Vector3D> r_1(a_search.size());
    std::vector<chemfiles::Vector3D> r_2(a_search.size());

    for (size_t i = 0; i < a_search.size(); ++i) {
        r_1[i] = a[a_search[i]];
        r_2[i] = b[a_native[i]];
    }

    int ier;
    std::vector<double> w(a_search.size(), 1.0);
    chemfiles::Vector3D t;
    chemfiles::Matrix3D u = chemfiles::Matrix3D::zero();
    lemon::tmalign::kabsch(w, r_1, r_2, a_search.size(), u, t, ier, 2.0, 2.0);
}
