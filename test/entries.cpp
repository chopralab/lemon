#include "lemon/entries.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <fstream>

TEST_CASE("Read entry file") {
    std::ifstream input_file("files/entries_10.idx");

    // Skip two lines
    std::string junk;
    std::getline(input_file, junk);
    std::getline(input_file, junk);

    auto vec = lemon::read_entry_file(input_file);
    
    CHECK(vec.size() == 10);
    CHECK(vec.count("100D") != 0);
    CHECK(vec.count("101D") != 0);
    CHECK(vec.count("101M") != 0);
    CHECK(vec.count("102D") != 0);
    CHECK(vec.count("102L") != 0);
    CHECK(vec.count("102M") != 0);
    CHECK(vec.count("103D") != 0);
    CHECK(vec.count("103L") != 0);
    CHECK(vec.count("103M") != 0);
    CHECK(vec.count("104D") != 0);
}

TEST_CASE("Read results file") {
    std::ifstream input_file("files/sam_results.idx");

    lemon::Entries vec;
    std::unordered_map<std::string, lemon::ResidueNameSet> rnms;
    lemon::read_entry_file(input_file, vec, rnms);
    CHECK(vec.size() == 2);
    CHECK(rnms.size() == 2);
    CHECK(rnms["4YND"].size() == 1);
    CHECK(rnms["5A1I"].size() == 2);
    CHECK(rnms["4YND"].count("4GQ") == 1);
    CHECK(rnms["5A1I"].count("ADN") == 1);
    CHECK(rnms["5A1I"].count("PPK") == 1);
}
