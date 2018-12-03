#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "lemon/hadoop.hpp"

#include <fstream>
#include <mutex>

#include "lemon/count.hpp"
#include "lemon/parallel.hpp"

TEST_CASE("Read single MMTF Sequence File") {
    std::ifstream hadoop_file("files/rcsb_hadoop/hadoop", std::istream::binary);
    lemon::Hadoop sequence(hadoop_file);
    auto result = sequence.next();

    // We've only got the one MMTF file
    CHECK(!sequence.has_next());

    CHECK(std::string(result.first.data() + 1, 4) == "1DZE");

    chemfiles::Trajectory traj(std::move(result.second), "MMTF/GZ");
    auto frame = traj.read();
}

TEST_CASE("Read multiple MMTF Sequence File") {
    std::ifstream hadoop_file("files/rcsb_hadoop/hadoop_multiple",
                              std::istream::binary);
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

TEST_CASE("Use Hadoop Run - no collector") {
    boost::filesystem::path p("files/rcsb_hadoop");
    std::map<std::string, size_t> counts;
    std::mutex test_mutex;

    auto worker = [&counts, &test_mutex](const chemfiles::Frame& complex,
                                         const std::string& pdbid) {
        auto result = lemon::count::bioassemblies(complex);
        std::lock_guard<std::mutex> guard(test_mutex);
        counts[pdbid] = result;
    };

    lemon::run_parallel(worker, p, 2);
    CHECK(counts.size() == 5);
}

TEST_CASE("Use Hadoop Run - with collector") {
    boost::filesystem::path p("files/rcsb_hadoop");

    auto worker = [](const chemfiles::Frame& complex,
                     const std::string&) {
        lemon::ResidueNameCount resn_counts;
        lemon::count::residues(complex, resn_counts);
        return resn_counts;
    };

    lemon::map_combine<lemon::ResidueNameCount> combiner;
    lemon::ResidueNameCount collector;

    lemon::run_parallel(worker, combiner, p, collector, 2);
    CHECK(collector.size() == 36);
}

TEST_CASE("Provide an invalid directory to Hadoop run") {
    CHECK_THROWS_AS(lemon::read_hadoop_dir({"/nodir/"}), std::runtime_error);
    CHECK_THROWS_AS(lemon::read_hadoop_dir({"."}), std::runtime_error);
    CHECK_THROWS_AS(lemon::read_hadoop_dir({"files/entry_10/1/0"}), std::runtime_error);
}
