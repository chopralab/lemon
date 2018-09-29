#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <fstream>

#include "lemon/count.hpp"
#include "lemon/hadoop.hpp"

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

    auto worker = [&counts](const chemfiles::Frame& complex,
                            const std::string& pdbid) {
        auto result = lemon::count_bioassemblies(complex);
        counts[pdbid] = result;
    };
    lemon::run_hadoop(worker, p);
    CHECK(counts.size() == 5);
}

TEST_CASE("Use Hadoop Run - with collector") {
    boost::filesystem::path p("files/rcsb_hadoop");
    std::map<std::string, size_t> counts;

    auto worker = [&counts](const chemfiles::Frame& complex,
                            const std::string&) {
        lemon::ResidueNameCount resn_counts;
        lemon::count_residues(complex, resn_counts);
        return resn_counts;
    };

    lemon::ResidueNameCount collector;

    lemon::run_hadoop(worker, p, collector);
    std::cout << "DERP{" << std::endl;
    CHECK(collector.size() == 36);
}

TEST_CASE("Provide an invalid directory to Hadoop run") {
    auto worker = [](const chemfiles::Frame&, const std::string&) {};

    CHECK_THROWS_AS(lemon::run_hadoop(worker, {"."}), std::runtime_error);
    CHECK_THROWS_AS(lemon::run_hadoop(worker, {"files/entry_10/1/0"}),
                    std::runtime_error);
}
