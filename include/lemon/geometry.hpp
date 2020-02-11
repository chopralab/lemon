#ifndef LEMON_GEOMETRY_HPP
#define LEMON_GEOMETRY_HPP

#include "lemon/constants.hpp"
#include "lemon/residue_name.hpp"
#include "lemon/external/gaurd.hpp"

#include <string>

LEMON_EXTERNAL_FILE_PUSH
#include "chemfiles/Frame.hpp"
LEMON_EXTERNAL_FILE_POP

namespace lemon {

//! Functions to retreive names of various bonds, angles, and dihedrals
namespace geometry {

//! Error to handle odd geometries
class geometry_error : public std::logic_error {
  public:
    explicit geometry_error(const std::string& what_arg)
        : std::logic_error(what_arg) {}
    explicit geometry_error(const char* what_arg)
        : std::logic_error(what_arg) {}
};

//! Is special residue
inline std::string is_special_residue(const std::string& resn,
                                      const ResidueNameSet& resn_names) {
    if (resn_names.count(resn) != 0) {
        return resn + '_';
    }

    return "";
}

//! For obtaining information about peptide geometry.
namespace protein {

//! Obtain the name of a `bond` present in a peptide residue
//!
//! This function returns the 'name' of a given `bond` in the form
//! {residue}_{atom1}_{atom2} for bonds within a single residue not part of the
//! amino acid linker. For bonds within the amino acid linker, only the atom
//! names are returned, ('CA_C', 'CA_N', 'C_O'), unless the residue name is in
//! `special` (C_O not included)). The default special residues are proline
//! derivatives. For peptide bonds and disulphide bridges, 'peptide_bond' and
//! 'SG_SG' are returned. Other types of inter-residue bonds result in a
//! geometry error being thrown. \param [in] entry Structure containing the
//! peptide residues of interest \param [in] bond Bond which name is going to be
//! obtained \param [in] special Residue names that should be handled
//! differently \return The name of the bond.
inline std::string bond_name(const chemfiles::Frame& entry,
                             const chemfiles::Bond& bond,
                             const ResidueNameSet& special = proline_res) {
    const auto& atom1 = entry[bond[0]];
    const auto& atom2 = entry[bond[1]];

    if (atom1.type() == "H") {
        return atom2.name() + "_H";
    }

    if (atom2.type() == "H") {
        return atom1.name() + "_H";
    }

    const std::string latom =
        atom1.name() < atom2.name() ? atom1.name() : atom2.name();

    const std::string hatom =
        atom1.name() > atom2.name() ? atom1.name() : atom2.name();

    // The carbonyl bond should all be the same
    if (latom == "C" && hatom == "O") {
        return "C_O";
    }

    const auto& residue1 = entry.topology().residue_for_atom(bond[0]);
    const auto& residue2 = entry.topology().residue_for_atom(bond[1]);

    if (residue1 != residue2) {
        if (latom == "C" && hatom == "N") {
            return "peptide_bond";
        }

        if (latom == "SG" && hatom == "SG") {
            return "SG_SG";
        }

        // clang-format off
        throw geometry_error( "Unhandled inter-residue bond: " +
                              residue1->name() + "_" + atom1.name() + " " +
                              residue2->name() + "_" + atom2.name());
        // clang-format on
    }

    // Let's keep the aminoacid group all the same (except proline)
    std::string name = is_special_residue(residue1->name(), special);

    if ((latom == "CA" || hatom == "CA") && !(latom == "CB" || hatom == "CB")) {
        return name + latom + "_" + hatom;
    }

    return residue1->name() + "_" + latom + "_" + hatom;
}

//! Obtain the name of an angle present in a peptide residue
//!
//! This function returns the 'name' of a given angle in the form
//! {residue}_{atom1}_{atom2}_{atom3} for angles within a single residue not
//! part of the amino acid linker, namely: CA_C_O, N_C_O, and O_C_OXT.
//! Inter-residue angles not part of the linker are handled similarly, both with
//! the exception of proline derivatives. Other types of inter-residue angles
//! result in a geometry error being thrown. \param [in] entry Structure
//! containing the peptide residues of interest \param [in] angle angle which
//! name is going to be obtained \return The name of the angle.
inline std::string angle_name(const chemfiles::Frame& entry,
                              const chemfiles::Angle& angle) {
    const auto& atom1 = entry[angle[0]];
    const auto& atom2 = entry[angle[1]];
    const auto& atom3 = entry[angle[2]];

    if (atom1.type() == "H") {
        return atom3.name() + "_" + atom2.name() + "_H";
    }

    if (atom3.type() == "H") {
        return atom1.name() + "_" + atom2.name() + "_H";
    }

    const std::string& latom =
        atom1.name() < atom3.name() ? atom1.name() : atom3.name();

    const std::string& hatom =
        atom1.name() > atom3.name() ? atom1.name() : atom3.name();

    const std::string& catom = atom2.name();

    // The carbonyl angles should all be the same
    if (catom == "C") {
        // Does it involve the carbonyl?
        if (latom == "O" || hatom == "O") {
            if (latom == "CA" || hatom == "CA") {
                return "CA_C_O";
            }
            if (latom == "N" || hatom == "N") {
                return "N_C_O"; // Inter-residue
            }
            if (latom == "OXT" || hatom == "OXT") {
                return "O_C_OXT";
            }
            throw geometry_error("Odd angle centered on 'C' atom " + latom +
                                 "_" + hatom);
        }

        // Maybe proline logic should be added here?
        // Potentially inter-residue (CA_C_N) or intra(CA_C_OXT)
        if (hatom == "N" || hatom == "OXT") {
            return latom + "_" + catom + "_" + hatom;
        }
        throw geometry_error("Odd angle centered on 'C' atom " + latom + "_" +
                             hatom);
    }

    if ((latom == "C" || hatom == "C") && catom == "CA") {
        // Special handling of the intra aminoacid group.
        // Another possibility is C_CA_CB, which is residue specific
        if ((hatom == "N") || (latom == "N")) {
            return "C_CA_N";
        }
    }

    if ((latom == "C" || hatom == "C") && catom == "N") {
        if (latom == "CA" || hatom == "CA") {
            return "C_N_CA";
        }
        // Proline-like!
        if (latom == "CD" || hatom == "CD") {
            return entry.topology().residue_for_atom(angle[1])->name() +
                   "_C_N_CD";
        }
        throw geometry_error("Unhandled common residue angle: " + latom + "_" +
                             catom + "_" + hatom);
    }

    // Check for sulfide bridges
    if (catom == "SG" && latom == "CB" && hatom == "SG") {
        return "CB_SG_SG";
    }

    const auto& residue1 = entry.topology().residue_for_atom(angle[0]);
    const auto& residue3 = entry.topology().residue_for_atom(angle[2]);

    // clang-format off
    if (residue1 != residue3) {
        const auto& residue2 = entry.topology().residue_for_atom(angle[1]);
        throw geometry_error("Unhandled inter-residue angle: " +
                             residue1->name() + "_" + atom1.name() + " " +
                             residue2->name() + "_" + atom2.name() + " " +
                             residue3->name() + "_" + atom3.name());
    }
    // clang-format on

    return residue1->name() + "_" + latom + "_" + catom + "_" + hatom;
}

//! Obtain the name of a dihedral present in a peptide residue
//!
//! This function returns the 'name' of a given dihedral in the form
//! {residue}_{atom1}_{atom2}_{atom3}_{atom4} for dihedrals within a single
//! residue not part of the amino acid linker. Namely: CA_C_N_CA, N_C_CA_N,
//! C_CA_N_C, N_CA_C_O, and CA_N_C_O. These are tagged with the residue name
//! if the residue is a proline derivative.
//! \param [in] entry Structure containing the peptide residues of interest
//! \param [in] dihedral Dihedral which name is going to be obtained
//! \param [in] special Residue names that should be handled differently
//! \return The name of the dihedral.
inline std::string dihedral_name(const chemfiles::Frame& entry,
                                 const chemfiles::Dihedral& dihedral,
                                 const ResidueNameSet& special = proline_res) {
    const auto& atom1 = entry[dihedral[0]];
    const auto& atom2 = entry[dihedral[1]];
    const auto& atom3 = entry[dihedral[2]];
    const auto& atom4 = entry[dihedral[3]];

    std::string llatom; // Left most atom
    std::string lcatom; // Center left atom
    std::string hcatom; // Center right
    std::string hhatom; // Right most

    if (atom1.name() < atom4.name()) {
        llatom = atom1.name();
        lcatom = atom2.name();
        hcatom = atom3.name();
        hhatom = atom4.name();
    } else if (atom4.name() < atom1.name()) {
        llatom = atom4.name();
        lcatom = atom3.name();
        hcatom = atom2.name();
        hhatom = atom1.name();
    } else { // They are equal
        llatom = hhatom = atom1.name();
        if (atom2.name() < atom3.name()) {
            lcatom = atom2.name();
            hcatom = atom3.name();
        } else {
            lcatom = atom3.name();
            hcatom = atom2.name();
        }
    }

    // We may need to be name specific instead of type specific here
    if (atom1.type() == "H") {
        return atom1.name() + "_" + atom2.name() + "_" + atom3.name() + "_H";
    }

    if (atom4.type() == "H") {
        return atom4.name() + "_" + atom3.name() + "_" + atom2.name() + "_H";
    }

    const auto& cresidue = entry.topology().residue_for_atom(dihedral[1]);

    std::string name = is_special_residue(cresidue->name(), special);

    // Peptide-bond, alpha-carbon to alpha carbon
    if (llatom == "CA" && hhatom == "CA") {
        return name + "CA_C_N_CA";
    }

    // Inter-residue dihedral across the entire residue
    // Nitrogen to Nitrogen
    if (llatom == "N" && hhatom == "N") {
        return name + "N_C_CA_N";
    }

    // Inter-residue dihedral across the entire residue
    // Carbonyl-carbon to carbonyl carbon
    if (llatom == "C" && hhatom == "C") {
        return name + "C_CA_N_C";
    }

    // Inter-residue dihedral for the peptide to N-CA bond
    // The version of this to the N-H bond is covered previously
    // The case for alpha-carbon to alpha-carbon is already covered
    if ((llatom == "N" && hhatom == "O") || (llatom == "O" && hhatom == "N")) {
        return name + "N_CA_C_O";
    }

    // Inter-residue dihedral for the peptide bond
    // The version of this to the N-H bond is covered previously
    if ((lcatom == "C" && hcatom == "N") || (lcatom == "N" && hcatom == "C")) {
        return name + "CA_N_C_O";
    }

    // Disulphide
    if (lcatom == "SG" && hcatom == "SG") {
        return "CB_SG_SG_CB";
    }

    if (hhatom == "SG" && hcatom == "SG") {
        return "CA_CB_SG_SG";
    }

    return cresidue->name() + "_" + llatom + "_" + lcatom + "_" + hcatom + "_" +
           hhatom;
}

//! Obtain the name of an improper dihedral present in a peptide residue
//!
//! This function returns the 'name' of a given improper in the form
//! {residue}_{atom1}_{atom2}_{atom3}_{atom4} for impropers within a single
//! residue not part of the amino acid linker (CA_C_N_O).
//! \param [in] entry Structure containing the peptide residues of interest
//! \param [in] improper Improper dihedral which name is going to be obtained
//! \param [in] special Residue names that should be handled differently
//! \return The name of the improper.
inline std::string improper_name(const chemfiles::Frame& entry,
                                 const chemfiles::Improper& improper,
                                 const ResidueNameSet& special = proline_res) {
    const auto& atom1 = entry[improper[0]];
    const auto& atom2 = entry[improper[1]];
    const auto& atom3 = entry[improper[2]];
    const auto& atom4 = entry[improper[3]];

    auto min_str = [](const std::string& s1, const std::string& s2) {
        return (s1 < s2) ? s1 : s2;
    };

    auto max_str = [](const std::string& s1, const std::string& s2) {
        return (s1 > s2) ? s1 : s2;
    };

    const std::string& latom =
        min_str(atom1.name(), min_str(atom3.name(), atom4.name()));

    const std::string& matom =
        max_str(min_str(atom1.name(), atom3.name()),
                min_str(max_str(atom1.name(), atom3.name()), atom4.name()));

    const std::string& hatom =
        max_str(atom1.name(), max_str(atom3.name(), atom4.name()));

    const std::string& catom = atom2.name();

    if (atom1.type() == "H") {
        return atom3.name() + "_" + atom2.name() + "_" + atom4.name() + "_H";
    }

    if (atom3.type() == "H") {
        return atom1.name() + "_" + atom2.name() + "_" + atom4.name() + "_H";
    }

    if (atom4.type() == "H") {
        return atom1.name() + "_" + atom2.name() + "_" + atom3.name() + "_H";
    }

    // The carbonyl angles should all be the same
    if (catom == "C") {
        // It must be the peptide bond!
        // as CA_C_N_H is already handled
        return "CA_C_N_O";
    }

    const auto& cresidue = entry.topology().residue_for_atom(improper[1]);

    std::string name = is_special_residue(cresidue->name(), special);

    if (catom == "N") {
        return name + latom + "_" + catom + "_" + matom + "_" + hatom;
    }

    return cresidue->name() + "_" + latom + "_" + catom + "_" + matom + "_" +
           hatom;
}
} // namespace protein
} // namespace geometry
} // namespace lemon

#endif
