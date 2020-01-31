#include "lemon/tmalign.hpp"

#include "chemfiles/Topology.hpp"
#include "chemfiles/Trajectory.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using Catch::Detail::Approx;

TEST_CASE("Kabsch") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1YT9.mmtf.gz", 'r');
    auto frame2 = traj.read();

    std::vector<size_t> a_search;
    std::vector<size_t> a_native;

    lemon::tmalign::find_operlapping_residues(frame1, frame2,
                                              a_search, a_native,
                                              "", "");

    auto a = frame1.positions();
    auto b = frame2.positions();

    std::vector<chemfiles::Vector3D> r_1(a_search.size());
    std::vector<chemfiles::Vector3D> r_2(a_search.size());

    for (size_t i = 0; i < a_search.size(); ++i) {
        r_1[i] = a[a_search[i]];
        r_2[i] = b[a_native[i]];
    }

    auto affine = lemon::kabsch(r_1, r_2);

    lemon::align(b, affine);
}

TEST_CASE("tmalign::TMscore for 1OQ5") {
    auto traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score ==  Approx(1));
    CHECK(tm.aligned == 256);
}

TEST_CASE("tmalign::TMscore for 1AAQ") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score ==  Approx(1));
    CHECK(tm.aligned == 199);
}

TEST_CASE("tmalign::TMscore for 1OQ5/1AAQ") {
    auto traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score == Approx(0.100864).epsilon(1e-3));
    CHECK(tm.aligned == 96);
}

TEST_CASE("tmalign::TMscore for 1AAQ/1OQ5") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score == Approx(0.0863468016).epsilon(1e-4));
    CHECK(tm.aligned == 96);
}

TEST_CASE("tmalign::TMscore for 1AAQ/4XUF") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score == Approx(0.0));
    CHECK(tm.aligned == 0);
}

TEST_CASE("tmalign::TMscore for 1AAQ/1YT9") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1YT9.mmtf.gz", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score == Approx(0.9762572135).epsilon(1e-4));
    CHECK(tm.aligned == 198);
}

TEST_CASE("tmalign::TMscore for 1K21/2V3O") {
    auto traj = chemfiles::Trajectory("files/1K21.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/2V3O.mmtf.gz", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score == Approx(0.9068108101).epsilon(1e-4));
    CHECK(tm.aligned == 276);
}

TEST_CASE("tmalign::TMscore for 1K21/1ZRB") {
    auto traj = chemfiles::Trajectory("files/1K21.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1ZRB.mmtf.gz", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score == Approx(0.8324166043).epsilon(1e-4));
    CHECK(tm.aligned == 248);
}

TEST_CASE("tmalign::TMscore for 1K21/1WAY") {
    auto traj = chemfiles::Trajectory("files/1K21.mmtf.gz", 'r');
    auto frame1 = traj.read();

    traj = chemfiles::Trajectory("files/1WAY.mmtf.gz", 'r');
    auto frame2 = traj.read();

    auto tm = lemon::tmalign::TMscore(frame1, frame2);
    CHECK(tm.score == Approx(0.8201441267).epsilon(1e-4));
    CHECK(tm.aligned == 248);
}
