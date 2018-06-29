#include "benchmarker/parse.hpp"
#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace benchmarker;

TEST_CASE("Residue Count") {
    ResnCount counter;
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
    //CHECK(counter["TRP"] == 0);
    CHECK(counter["MET"] == 4);

    CHECK(counter["HOH"] == 1);
    CHECK(counter["CA"]  == 6);
    CHECK(counter["K"]   == 1);
}