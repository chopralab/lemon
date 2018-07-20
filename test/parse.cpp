#include "lemon/parse.hpp"
#include <chemfiles.hpp>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace lemon;

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

TEST_CASE("Assembly count") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto test = traj.read();

    CHECK(count_bioassemblies(test) == 1);

    traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    test = traj.read();

    CHECK(count_bioassemblies(test) == 2);

    traj = chemfiles::Trajectory("files/2WTL.mmtf.gz", 'r');
    test = traj.read();

    CHECK(count_bioassemblies(test) == 1);
}
