#include "benchmarker/parse.hpp"
#include <chemfiles.hpp>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace benchmarker;

TEST_CASE("Residue Name") {
    REQUIRE(sizeof(ResidueName) == sizeof(char) * 3);

    // String cannot be over 3 digits, or have no digits
    CHECK_THROWS_AS(ResidueName("ASDF"), std::length_error);
    CHECK_THROWS_AS(ResidueName(""), std::length_error);

    // Must be alpha numeric
    CHECK_THROWS_AS(ResidueName("+"), std::range_error);
    CHECK_THROWS_AS(ResidueName(" "), std::range_error);

    CHECK(ResidueName("ABC") == std::string("ABC"));
    CHECK(ResidueName("AB") == std::string("AB"));
    CHECK(ResidueName("A") == std::string("A"));

    CHECK(std::string("A") == ResidueName("A"));

    CHECK(ResidueName("ABD") != std::string("ABC"));
    CHECK(ResidueName("AB") != std::string("ABC"));
    CHECK(ResidueName("A") != std::string("AB"));

    CHECK(std::string("AB") != ResidueName("ABC"));
    CHECK(std::string("A") != ResidueName("AB"));

    CHECK(ResidueName("000").hash() == 0);
    CHECK(ResidueName("100").hash() == 1);
    CHECK(ResidueName("900").hash() == 9);
    CHECK(ResidueName("A00").hash() == 10);
    CHECK(ResidueName("Z00").hash() == 35);
    // Note: a hash of 36 is impossible as ' 00' is not valid

    CHECK(ResidueName("010").hash() == 37);
    CHECK(ResidueName("110").hash() == 38);
    CHECK(ResidueName("210").hash() == 39);
    CHECK(ResidueName("310").hash() == 40);

    CHECK(ResidueName("020").hash() == 74);
    CHECK(ResidueName("120").hash() == 75);

    CHECK(ResidueName("0Z0").hash() == 37 * 35);
    CHECK(ResidueName("ZZ0").hash() == 37 * 35 + 35);
    CHECK(ResidueName("ZZZ").hash() == 35 + 35 * 37 + 35 * 37 * 37);

    CHECK(ResidueName("01").hash() == 0 + 1 * 37 + 36 * 37 * 37);
    CHECK(ResidueName("11").hash() == 1 + 1 * 37 + 36 * 37 * 37);
    CHECK(ResidueName("02").hash() == 0 + 2 * 37 + 36 * 37 * 37);

    CHECK(ResidueName("0").hash() == 0 + 36 * 37 + 36 * 37 * 37);
    CHECK(ResidueName("1").hash() == 1 + 36 * 37 + 36 * 37 * 37);
    CHECK(ResidueName("Z").hash() == 35 + 36 * 37 + 36 * 37 * 37);
}

TEST_CASE("Residue Count") {
    ResidueNameCount counter;
    auto traj = chemfiles::Trajectory("files/5w1d.pdb", 'r');
    auto test1 = traj.read();


    retreive_residue_counts(test1, counter);

    CHECK(counter.size() == 22);

    CHECK(counter["ALA"] == 28);
    CHECK(counter["ARG"] == 19);
    CHECK(counter["ASN"] == 24);
    CHECK(counter["ASP"] == 27);
    CHECK(counter["CYS"] == 1);
    CHECK(counter["GLN"] == 14);
    CHECK(counter["GLU"] == 25);
    CHECK(counter["GLY"] == 26);
    CHECK(counter["HIS"] == 6);
    CHECK(counter["ILE"] == 28);
    CHECK(counter["MET"] == 4);
    CHECK(counter["LEU"] == 37);
    CHECK(counter["LYS"] == 5);
    CHECK(counter["PHE"] == 12);
    CHECK(counter["PRO"] == 25);
    CHECK(counter["SER"] == 27);
    CHECK(counter["THR"] == 37);
    CHECK(counter["TYR"] == 19);
    CHECK(counter["TRP"] == 0);
    CHECK(counter["MET"] == 4);

    CHECK(counter["HOH"] == 1);
    CHECK(counter["CA"]  == 6);
    CHECK(counter["K"]   == 1);

    traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto test2 = traj.read();

    retreive_residue_counts(test2, counter);
    CHECK(counter.size() == 24);

    // Check the simple updates
    CHECK(counter["HOH"] == 2);
    CHECK(counter["PSI"] == 1);
    CHECK(counter["TRP"] == 4);

    // No change
    CHECK(counter["CA"]  == 6);
    CHECK(counter["K"]   == 1);

    ResidueNameCount counter2;
    retreive_residue_counts(test1, counter2);
    retreive_residue_counts(test2, counter2);

    CHECK(counter2 == counter);
}

TEST_CASE("Residue to text") {
    auto res_name = ResidueName("ALA");
    std::stringstream ss;

    ss << res_name;
    CHECK(ss.str() == "ALA");

    std::stringstream ss2;
    ResidueNameCount rnc;
    rnc[res_name] = 75;
    ss2 << rnc;
    CHECK(ss2.str() == "ALA\t75\n");
}
