#include "lemon/residue_name.hpp"
#include <chemfiles.hpp>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace lemon;

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

TEST_CASE("Residue set") {
    ResidueNameSet rns {"ALA", "GLY", "ALA"};
    CHECK(rns.size() == 2);

    rns.insert("GLY");
    CHECK(rns.size() == 2);

    rns.insert("VAL");
    CHECK(rns.size() == 3);

    rns.insert("GLY");
    CHECK(rns.size() == 3);
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
