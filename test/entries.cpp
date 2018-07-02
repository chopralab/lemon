#include "benchmarker/entries.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <fstream>

TEST_CASE("Read entry file") {
    std::ifstream input_file("files/entries_10.idx");

    // Skip two lines
    std::string junk;
    std::getline(input_file, junk);
    std::getline(input_file, junk);

    std::vector<std::array<char, 4>> vec;
    vec.reserve(10);
    benchmarker::read_entry_file(input_file, vec);
    
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
