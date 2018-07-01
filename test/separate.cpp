#include "benchmarker/separate.hpp"
#include <chemfiles.hpp>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Separate") {
    auto traj = chemfiles::Trajectory("files/1AAQ.mmtf", 'r');
    auto frame = traj.read();

    chemfiles::Frame ligand;
    chemfiles::Frame protein;

    const auto& residues = frame.topology().residues();
    size_t i; 
    for (i = 0; i < residues.size(); ++i) {
        if (residues[i].name() == "PSI") {
            break;
        }
    }

    benchmarker::separate_protein_and_ligand(frame, residues[i], protein, ligand, 15);

    CHECK(ligand.size() == 41);
    CHECK(protein.size() == 998);
}
