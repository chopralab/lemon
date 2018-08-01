#include "lemon/hadoop.hpp"
#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <fstream>

TEST_CASE("Read single MMTF Sequence File") {
    std::ifstream hadoop_file("files/hadoop.seq");
    lemon::Hadoop sequence(hadoop_file);
    auto result = sequence.next();

    // We've only got the one MMTF file
    CHECK(!sequence.has_next());

    CHECK(std::string(result.first.data() + 1, 4) == "1DZE");

    chemfiles::Trajectory traj(std::move(result.second), "MMTF/GZ");
    auto frame = traj.read();
}

TEST_CASE("Read multiple MMTF Sequence File") {
    std::ifstream hadoop_file("files/hadoop_multiple.seq");
    lemon::Hadoop sequence(hadoop_file);
    size_t count = 0;
    while (sequence.has_next()) {
        auto result = sequence.next();
        chemfiles::Trajectory traj(std::move(result.second), "MMTF/GZ");
        auto frame = traj.read();
        ++count;
    }
    CHECK(count == 5);
}
