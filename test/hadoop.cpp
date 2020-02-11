#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "lemon/hadoop.hpp"

#include <fstream>
#include <mutex>

#include "lemon/count.hpp"
#include "lemon/parallel.hpp"
#include "lemon/launch.hpp"

TEST_CASE("Read single MMTF Sequence File") {
    std::ifstream hadoop_file("files/rcsb_hadoop/hadoop", std::istream::binary);
    lemon::Hadoop sequence(hadoop_file);
    auto result = sequence.next();

    // We've only got the one MMTF file
    CHECK(!sequence.has_next());

    CHECK(result.first == "1DZE");

    auto traj = chemfiles::Trajectory::memory_reader(
        result.second.data(),
        result.second.size(),
        "MMTF/GZ"
    );
    auto frame = traj.read();
}

TEST_CASE("Read multiple MMTF Sequence File") {
    std::ifstream hadoop_file("files/rcsb_hadoop/hadoop_multiple",
                              std::istream::binary);
    lemon::Hadoop sequence(hadoop_file);
    size_t count = 0;
    while (sequence.has_next()) {
        auto result = sequence.next();
        auto traj = chemfiles::Trajectory::memory_reader(
            result.second.data(),
            result.second.size(),
            "MMTF/GZ"
        );
        auto frame = traj.read();
        ++count;
    }
    CHECK(count == 5);
}

TEST_CASE("Use run_parallel") {
    std::string p("files/rcsb_hadoop");

    auto worker = [](const chemfiles::Frame& entry,
                     const std::string& /*unused*/) -> lemon::ResidueNameCount {
        lemon::ResidueNameCount resn_counts;
        lemon::count::residues(entry, resn_counts);
        return resn_counts;
    };

    lemon::ResidueNameCount totals;
    auto collector = lemon::map_combine<lemon::ResidueNameCount>(totals);

    lemon::run_parallel(worker, p, collector, 2);
    CHECK(totals.size() == 36);
}

TEST_CASE("Use run_parallel, but only for a entry") {
    std::string p("files/rcsb_hadoop");

    auto worker = [](const chemfiles::Frame& entry,
                     const std::string& /*unused*/) -> lemon::ResidueNameCount {
        lemon::ResidueNameCount resn_counts;
        lemon::count::residues(entry, resn_counts);
        return resn_counts;
    };

    lemon::ResidueNameCount totals;
    auto collector = lemon::map_combine<lemon::ResidueNameCount>(totals);

    std::unordered_set<std::string> e({"1DZE"});

    lemon::run_parallel(worker, p, collector, 2, e);
    CHECK(totals.size() == 27); // 1DZE has a lot of cofactors...
}

TEST_CASE("Use run_parallel, but skip a entry") {
    std::string p("files/rcsb_hadoop");

    auto worker = [](const chemfiles::Frame& entry,
                     const std::string& /*unused*/) -> lemon::ResidueNameCount {
        lemon::ResidueNameCount resn_counts;
        lemon::count::residues(entry, resn_counts);
        return resn_counts;
    };

    lemon::ResidueNameCount totals;
    auto collector = lemon::map_combine<lemon::ResidueNameCount>(totals);
    std::unordered_set<std::string> se({"1DZE"});

    lemon::run_parallel(worker, p, collector, 2, lemon::Entries(), se);
    CHECK(totals.size() == 28);
}

TEST_CASE("Provide an invalid directory to Hadoop run") {
    CHECK_THROWS_AS(lemon::read_hadoop_dir({"/nodir/"}), std::runtime_error&);
    CHECK_THROWS_AS(lemon::read_hadoop_dir({"."}), std::runtime_error&);
    CHECK_THROWS_AS(lemon::read_hadoop_dir({"files/entry_10/1/0"}), std::runtime_error&);
}
