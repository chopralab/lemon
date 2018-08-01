#include "lemon/hadoop.hpp"
#include <chemfiles.hpp>
#include "lemon/count.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <fstream>

TEST_CASE("Read single MMTF Sequence File") {
    std::ifstream hadoop_file("files/rcsb_hadoop/hadoop.seq");
    lemon::Hadoop sequence(hadoop_file);
    auto result = sequence.next();

    // We've only got the one MMTF file
    CHECK(!sequence.has_next());

    CHECK(std::string(result.first.data() + 1, 4) == "1DZE");

    chemfiles::Trajectory traj(std::move(result.second), "MMTF/GZ");
    auto frame = traj.read();
}

TEST_CASE("Read multiple MMTF Sequence File") {
    std::ifstream hadoop_file("files/rcsb_hadoop/hadoop_multiple.seq");
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

TEST_CASE("Use Hadoop Run") {
    boost::filesystem::path p("files/rcsb_hadoop");
    std::map<std::string, size_t> counts;

    auto worker = [&counts](const chemfiles::Frame& complex,
                              const std::string& pdbid) {

        auto result = lemon::count_bioassemblies(complex);
        counts[pdbid] = result;
        std::cout << "here" << std::endl;
    };
    lemon::run_hadoop(worker, p);

    CHECK(counts.size() == 5);
}
