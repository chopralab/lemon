#include "lemon/geometry.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/Trajectory.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

namespace lg = lemon::geometry;

TEST_CASE("Geometry for 1OQ5") {
    auto traj = chemfiles::Trajectory("files/1OQ5.mmtf.gz", 'r');
    auto frame1 = traj.read();

    // Peptide bond
    CHECK(lg::protein::bond_name(frame1, {2, 10}) == "peptide_bond");

    // Bonds part of the linker
    CHECK(lg::protein::bond_name(frame1, {0, 1}) == "CA_N");
    CHECK(lg::protein::bond_name(frame1, {1, 2}) == "C_CA");
    CHECK(lg::protein::bond_name(frame1, {2, 3}) == "C_O");

    // Not part of the linker
    CHECK(lg::protein::bond_name(frame1, {1, 4}) == "HIS_CA_CB");

    // Proline!
    CHECK(lg::protein::bond_name(frame1, {71, 72}) == "PRO_CA_N");
    CHECK(lg::protein::bond_name(frame1, {73, 72}) == "PRO_C_CA");

    // Mistake!
    CHECK_THROWS(lg::protein::bond_name(frame1, {25, 100}));

    CHECK(lg::protein::angle_name(frame1, {1, 2, 3}) == "CA_C_O");
    CHECK(lg::protein::angle_name(frame1, {1, 2, 10}) == "CA_C_N");
    CHECK(lg::protein::angle_name(frame1, {3, 2, 10}) == "N_C_O");
    CHECK_THROWS(lg::protein::angle_name(frame1, {1, 2, 107}));
    CHECK_THROWS(lg::protein::angle_name(frame1, {3, 2, 107}));

    CHECK(lg::protein::angle_name(frame1, {0, 1, 2}) == "C_CA_N");
    CHECK(lg::protein::angle_name(frame1, {0, 1, 4}) == "HIS_CB_CA_N");
    CHECK(lg::protein::angle_name(frame1, {4, 1, 2}) == "HIS_C_CA_CB");

    CHECK(lg::protein::angle_name(frame1, {2, 10, 11}) == "C_N_CA");
    CHECK(lg::protein::angle_name(frame1, {73, 71, 77}) == "PRO_C_N_CD");
    CHECK_THROWS(lg::protein::angle_name(frame1, {2, 10, 86}));

    CHECK(lg::protein::improper_name(frame1, {0, 2, 3, 10}) == "CA_C_N_O");
    CHECK(lg::protein::improper_name(frame1, {37, 38, 36, 39}) == "TYR_CE1_CZ_CE2_OH");
    CHECK(lg::protein::improper_name(frame1, {140, 149, 150, 155}) == "PRO_C_N_CA_CD");

    CHECK(lg::protein::dihedral_name(frame1, {1, 2, 10, 11}) == "CA_C_N_CA");
    CHECK(lg::protein::dihedral_name(frame1, {0, 1, 2, 10}) == "N_C_CA_N");
    CHECK(lg::protein::dihedral_name(frame1, {2, 10, 11, 12}) == "C_CA_N_C");
    CHECK(lg::protein::dihedral_name(frame1, {0, 1, 2, 3}) == "N_CA_C_O");
    CHECK(lg::protein::dihedral_name(frame1, {11, 10, 2, 3}) == "CA_N_C_O");
    CHECK(lg::protein::dihedral_name(frame1, {11, 14, 15, 16}) == "TRP_CA_CB_CG_CD1");
}

TEST_CASE("Geometry for 1D7D") {
    auto traj = chemfiles::Trajectory("files/1D7D.mmtf.gz", 'r');
    auto frame1 = traj.read();

    CHECK(lg::protein::bond_name(frame1, {917, 936}) == "SG_SG");
    CHECK(lg::protein::angle_name(frame1, {916, 917, 936}) == "CB_SG_SG");
    CHECK(lg::protein::angle_name(frame1, {917, 936, 935}) == "CB_SG_SG");
    CHECK(lg::protein::dihedral_name(frame1, {916, 917, 936, 935}) == "CB_SG_SG_CB");
    CHECK(lg::protein::dihedral_name(frame1, {915, 916, 917, 936}) == "CA_CB_SG_SG");
}
