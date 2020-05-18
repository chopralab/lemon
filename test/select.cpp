#include "lemon/select.hpp"
#include <chemfiles.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Select PSI from 1AAQ") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame = traj.read();

    const auto res = lemon::select::small_molecules(frame);
    CHECK(res.size() == 1);

    size_t id = *res.begin();
    CHECK(frame.topology().residues()[id].name() == "PSI");
}

TEST_CASE("Select P30 from 4XUF") {
    auto traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame = traj.read();

    const auto res = lemon::select::small_molecules(frame);
    CHECK(res.size() == 2);

    size_t id = *res.begin();
    CHECK(frame.topology().residues()[id].name() == "P30");
}

TEST_CASE("Select FE from 2WTL") {
    auto traj = chemfiles::Trajectory("files/2WTL.mmtf.gz", 'r');
    auto frame = traj.read();

    const auto res = lemon::select::metal_ions(frame);
    CHECK(res.size() == 12);
}

TEST_CASE("Select CD from 1D7D") {
    auto traj = chemfiles::Trajectory("files/1D7D.mmtf.gz", 'r');
    auto frame = traj.read();

    const auto res = lemon::select::metal_ions(frame);
    CHECK(res.size() == 6);

    std::unordered_map<std::string, size_t> metals;
    const auto& residues = frame.topology().residues();

    for (auto res_id : res) {
        auto atom_id = *(residues[res_id].begin());
        auto iter = metals.find(frame[atom_id].type());
        if (iter == metals.end()) {
            metals[frame[atom_id].type()] = 1;
            continue;
        }
        ++iter->second;
    }

    CHECK(metals["Cd"] == 6);
}

TEST_CASE("Select HOH from 1OQ5") {
    auto traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame = traj.read();

    const auto res = lemon::select::specific_residues(frame, {"HOH"});
    CHECK(res.size() == 233);

    // select by ID
    const auto res2 = lemon::select::residue_ids(frame, {701});
    const auto& residue_701 = frame.topology().residues()[*res2.begin()];
    CHECK(residue_701.size() == 26);
    CHECK(residue_701.name() == "CEL");
}

TEST_CASE("Select nothing from 2WTL") {
    auto traj = chemfiles::Trajectory("files/2WTL.mmtf.gz", 'r');
    auto frame = traj.read();

    auto res = lemon::select::small_molecules(frame);
    CHECK(res.size() == 3); //3 UNLs exist in the structure
}

TEST_CASE("Select peptides") {
    auto traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame = traj.read();

    auto peptides = lemon::select::peptides(frame);
    CHECK(peptides.size() == 544);

    traj = chemfiles::Trajectory("files/entry_10/1/0/100D.mmtf.gz", 'r');
    frame = traj.read();

    peptides = lemon::select::peptides(frame);
    CHECK(peptides.size() == 0);
}

TEST_CASE("Select the chain name") {
    auto traj = chemfiles::Trajectory("files/4XUF.mmtf.gz", 'r');
    auto frame = traj.read();

    auto chain_a = lemon::select::residue_property(frame, "chainname", chemfiles::Property("A"));

    CHECK(chain_a.size() == 273);
}
