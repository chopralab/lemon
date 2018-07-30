#include <iostream>

#include "lemon/entries.hpp"
#include "lemon/count.hpp"
#include "lemon/run.hpp"
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
    lemon::call_multithreaded(worker, vec, 1, 1);

    std::map<std::string, size_t> counts_1(counts);
    counts.clear();

    lemon::call_multithreaded(worker, vec, 2, 1);

    std::map<std::string, size_t> counts_2(counts);
    counts.clear();

    lemon::call_multithreaded(worker, vec, 2, 5);

    std::map<std::string, size_t> counts_3(counts);
    counts.clear();

    CHECK(counts_1 == counts_2);
    CHECK(counts_2 == counts_3);
}

TEST_CASE("Residue names") {
    auto worker1 = [](const chemfiles::Frame&,
                      const std::string&) {
        lemon::ResidueName("+");
    };

    auto worker2 = [](const chemfiles::Frame&,
                      const std::string&) {
        lemon::ResidueName("ABCD");
    };

    lemon::PDBIDVec vec;
    vec.reserve(10);
    lemon::read_entry_file("files/entries_10.idx", vec);

    // Exceptions should not escape
    lemon::call_multithreaded(worker1, vec, 1, 1);
    lemon::call_multithreaded(worker2, vec, 1, 1);
}

TEST_CASE("Other issues") {
    auto worker1 = [](const chemfiles::Frame&,
                      const std::string&) {
        throw std::logic_error("ASDF");
    };

    auto worker2 = [](const chemfiles::Frame&,
                      const std::string&) {
        throw 5;
    };

    lemon::PDBIDVec vec;
    vec.reserve(10);
    lemon::read_entry_file("files/entries_10.idx", vec);

    // Exceptions should not escape
    lemon::call_multithreaded(worker1, vec, 1, 1);
    lemon::call_multithreaded(worker2, vec, 1, 1);
}
