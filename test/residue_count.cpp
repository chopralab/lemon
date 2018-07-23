#include "lemon/residue_name.hpp"
#include <chemfiles.hpp>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Residue Name") {
    REQUIRE(sizeof(lemon::ResidueName) == sizeof(char) * 3);

    // String cannot be over 3 digits, or have no digits
    CHECK_THROWS_AS(lemon::ResidueName("ASDF"), std::length_error);
    CHECK_THROWS_AS(lemon::ResidueName(std::string("")), std::length_error);

    // Must be alpha numeric
    CHECK_THROWS_AS(lemon::ResidueName("+"), std::range_error);
    CHECK_THROWS_AS(lemon::ResidueName(std::string(" ")), std::range_error);

    CHECK(lemon::ResidueName("ABC") == std::string("ABC"));
    CHECK(lemon::ResidueName("AB") == std::string("AB"));
    CHECK(lemon::ResidueName("A") == std::string("A"));

    CHECK(std::string("A") == lemon::ResidueName("A"));

    CHECK(lemon::ResidueName("ABD") != std::string("ABDC"));
    CHECK(lemon::ResidueName("ABD") != std::string("ABC"));
    CHECK(lemon::ResidueName("AB") != std::string("ABC"));
    CHECK(lemon::ResidueName("A") != std::string("AB"));

    CHECK(std::string("AB") != lemon::ResidueName("ABC"));
    CHECK(std::string("A") != lemon::ResidueName("AB"));

    CHECK(lemon::ResidueName("000").hash() == 0);
    CHECK(lemon::ResidueName("100").hash() == 1);
    CHECK(lemon::ResidueName("900").hash() == 9);
    CHECK(lemon::ResidueName("A00").hash() == 10);
    CHECK(lemon::ResidueName("Z00").hash() == 35);
    // Note: a hash of 36 is impossible as ' 00' is not valid

    CHECK(lemon::ResidueName("010").hash() == 37);
    CHECK(lemon::ResidueName("110").hash() == 38);
    CHECK(lemon::ResidueName("210").hash() == 39);
    CHECK(lemon::ResidueName("310").hash() == 40);

    CHECK(lemon::ResidueName("020").hash() == 74);
    CHECK(lemon::ResidueName("120").hash() == 75);

    CHECK(lemon::ResidueName("0Z0").hash() == 37 * 35);
    CHECK(lemon::ResidueName("ZZ0").hash() == 37 * 35 + 35);
    CHECK(lemon::ResidueName("ZZZ").hash() == 35 + 35 * 37 + 35 * 37 * 37);

    CHECK(lemon::ResidueName("01").hash() == 0 + 1 * 37 + 36 * 37 * 37);
    CHECK(lemon::ResidueName("11").hash() == 1 + 1 * 37 + 36 * 37 * 37);
    CHECK(lemon::ResidueName("02").hash() == 0 + 2 * 37 + 36 * 37 * 37);

    CHECK(lemon::ResidueName("0").hash() == 0 + 36 * 37 + 36 * 37 * 37);
    CHECK(lemon::ResidueName("1").hash() == 1 + 36 * 37 + 36 * 37 * 37);
    CHECK(lemon::ResidueName("Z").hash() == 35 + 36 * 37 + 36 * 37 * 37);
}

TEST_CASE("Residue Name Counts") {
    lemon::ResidueNameCount rnc;
    rnc["ASP"] = 1;
    rnc["TRP"] = 2;
    CHECK(rnc["ASP"] == 1);
    CHECK(rnc["TRP"] == 2);

    lemon::ResidueNameCount rnc_2;
    rnc_2["ASP"] = 1;
    rnc_2["LYS"] = 2;
    rnc += rnc_2;

    CHECK(rnc["ASP"] == 2);
    CHECK(rnc["TRP"] == 2);
    CHECK(rnc["LYS"] == 2);
}

TEST_CASE("Residue set") {
    lemon::ResidueNameSet rns {"ALA", "GLY", "ALA"};
    CHECK(rns.size() == 2);

    rns.insert("GLY");
    CHECK(rns.size() == 2);

    rns.insert("VAL");
    CHECK(rns.size() == 3);

    rns.insert("GLY");
    CHECK(rns.size() == 3);
}

TEST_CASE("Residue to text") {
    auto res_name = lemon::ResidueName("ALA");
    std::stringstream ss;

    ss << res_name;
    CHECK(ss.str() == "ALA");

    std::stringstream ss2;
    lemon::ResidueNameCount rnc;
    rnc[res_name] = 75;
    ss2 << rnc;
    CHECK(ss2.str() == "\tALA\t75");
}
