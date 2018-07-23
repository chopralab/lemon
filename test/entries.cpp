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

    lemon::PDBIDVec vec;
    vec.reserve(10);
    lemon::read_entry_file(input_file, vec);
    
    CHECK(vec.size() == 10);
    CHECK(std::string(vec[0].data(), 4) == "100D");
    CHECK(std::string(vec[1].data(), 4) == "101D");
    CHECK(std::string(vec[2].data(), 4) == "101M");
    CHECK(std::string(vec[3].data(), 4) == "102D");
    CHECK(std::string(vec[4].data(), 4) == "102L");
    CHECK(std::string(vec[5].data(), 4) == "102M");
    CHECK(std::string(vec[6].data(), 4) == "103D");
    CHECK(std::string(vec[7].data(), 4) == "103L");
    CHECK(std::string(vec[8].data(), 4) == "103M");
    CHECK(std::string(vec[9].data(), 4) == "104D");
}

TEST_CASE("Read results file") {
    std::ifstream input_file("files/sam_results.idx");

    lemon::PDBIDVec vec;
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
