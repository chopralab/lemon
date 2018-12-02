#include "lemon/count.hpp"
#include <chemfiles.hpp>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Residue Count") {
    lemon::ResidueNameCount counter;
    auto traj = chemfiles::Trajectory("files/5w1d.pdb", 'r');
    auto test1 = traj.read();

    lemon::count::residues(test1, counter);

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

    lemon::count::residues(test2, counter);
    CHECK(counter.size() == 24);

    // Check the simple updates
    CHECK(counter["HOH"] == 2);
    CHECK(counter["PSI"] == 1);
    CHECK(counter["TRP"] == 4);

    // No change
    CHECK(counter["CA"]  == 6);
    CHECK(counter["K"]   == 1);

    lemon::ResidueNameCount counter2;
    lemon::count::residues(test1, counter2);
    lemon::count::residues(test2, counter2);

    CHECK(counter2 == counter);
}

TEST_CASE("Selective residue count") {
    lemon::ResidueNameCount counter;
    auto traj = chemfiles::Trajectory("files/5w1d.pdb", 'r');
    auto test1 = traj.read();

    lemon::count::residues(test1, {1, 2, 3, 4, 5, 6}, counter);

    auto s1 = lemon::count::print_residue_name_counts("JUNK", test1, {});
    CHECK(s1 == "");

    auto s2 = lemon::count::print_residue_name_counts("5w1d", test1, {2, 5});
    CHECK(s2== "5w1d\tTYR\t2\n");
}

TEST_CASE("Alternative Location count") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto test = traj.read();
    CHECK(lemon::count::altloc(test) == 0);

    traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    test = traj.read();
    CHECK(lemon::count::altloc(test) == 0);

    traj = chemfiles::Trajectory("files/2WTL.mmtf.gz", 'r');
    test = traj.read();
    CHECK(lemon::count::altloc(test) == 2);

    traj = chemfiles::Trajectory("files/entry_10/5/A/5A1I.mmtf.gz", 'r');
    test = traj.read();
    CHECK(lemon::count::altloc(test) == 4);
}

TEST_CASE("Assembly count") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto test = traj.read();
    CHECK(lemon::count::bioassemblies(test) == 1);

    traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    test = traj.read();
    CHECK(lemon::count::bioassemblies(test) == 2);

    traj = chemfiles::Trajectory("files/2WTL.mmtf.gz", 'r');
    test = traj.read();
    CHECK(lemon::count::bioassemblies(test) == 1);
}
