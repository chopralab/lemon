#include <iostream>

#include "lemon/entries.hpp"
#include "lemon/count.hpp"
#include "lemon/archive_run.hpp"
#include "lemon/select.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <boost/filesystem.hpp>

#include <fstream>

using namespace boost::filesystem;

TEST_CASE("Ensure threading works") {

    lemon::PDBIDVec vec;
    vec.reserve(10);
    lemon::read_entry_file("files/entries_10.idx", vec);

    std::map<std::string, size_t> counts;

    auto worker = [&counts](const chemfiles::Frame& complex,
                              const std::string& pdbid) {

        auto result = lemon::count_bioassemblies(complex);
        counts[pdbid] = result;
    };

    current_path("files/entry_10");
    lemon::run_archive(worker, vec, 1, 1);

    std::map<std::string, size_t> counts_1(counts);
    counts.clear();

    lemon::run_archive(worker, vec, 2, 1);

    std::map<std::string, size_t> counts_2(counts);
    counts.clear();

    lemon::run_archive(worker, vec, 2, 5);

    std::map<std::string, size_t> counts_3(counts);
    counts.clear();

    CHECK(counts_1 == counts_2);
    CHECK(counts_2 == counts_3);
}

