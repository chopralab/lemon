#include "lemon/separate.hpp"
#include "lemon/select.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/Trajectory.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Separate a ligand") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame = traj.read();

    chemfiles::Frame ligand;
    chemfiles::Frame protein;

    const auto psi_residue = lemon::select::specific_residues(frame, {"PSI"});

    // NOLINTNEXTLINE allow 15 to be magic
    lemon::separate::protein_and_ligand(frame, *psi_residue.begin(), 15, protein, ligand);

    CHECK(ligand.size() == 41);
    CHECK(protein.size() == 998);

    CHECK(protein.topology().residues()[0].get("is_standard_pdb")->as_bool());
}

TEST_CASE("Separate ligands") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame = traj.read();

    chemfiles::Frame ligand;
    chemfiles::Frame protein;

    const auto psi_residue = lemon::select::specific_residues(frame, {"PSI", "HOH"});

    // NOLINTNEXTLINE allow 15 to be magic
    lemon::separate::protein_and_ligands(frame, psi_residue, 15, protein, ligand);

    CHECK(ligand.size() == 42);
    CHECK(protein.size() == 997);

    CHECK(protein.topology().residues()[0].get("is_standard_pdb")->as_bool());
}
