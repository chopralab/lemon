#include <iostream>
#include <sstream>
#include <unordered_map>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/count.hpp"
#include "lemon/entries.hpp"
#include "lemon/hadoop.hpp"
#include "lemon/options.hpp"
#include "lemon/prune.hpp"
#include "lemon/select.hpp"
#include "lemon/separate.hpp"

// typedefs for binned data
typedef std::pair<std::string, int> BondDihedralBin;
typedef std::map<BondDihedralBin, size_t> DihedralCounts;

inline std::string get_dihedral_name(const chemfiles::Frame& complex,
                                     const chemfiles::Dihedral& dihedral) {
    const auto& atom1 = complex[dihedral[0]];
    const auto& atom2 = complex[dihedral[1]];
    const auto& atom3 = complex[dihedral[2]];
    const auto& atom4 = complex[dihedral[3]];

    std::string llatom, lcatom, hcatom, hhatom;
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
    } else {
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
        return atom1.type() + "_" + atom2.type() + "_" + atom3.type() + "_H";
    }

    if (atom4.type() == "H") {
        return atom4.type() + "_" + atom3.type() + "_" + atom2.type() + "_H";
    }

    const auto& cresidue = complex.topology().residue_for_atom(dihedral[1]);

    std::string name;
    // Proline is ... special
    if (cresidue->name() == "PRO") {
        name = "PRO_";
    }

    // Same goes for hydroxyproline
    if (cresidue->name() == "HYP") {
        name = "HYP_";
    }

    // and pyroglutamic acid
    if (cresidue->name() == "PCA") {
        name = "PCA_";
    }

    // Peptide-bond, alpha-carbon to alpha carbon
    if (llatom == "CA" && hhatom == "CA") {
        return "CA_C_N_CA";
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
    if ((llatom == "N" && hhatom == "O") ||
        (llatom == "O" && hhatom == "N")) {
        return name + "N_CA_C_O";
    }

    // Inter-residue dihedral for the peptide bond
    // The version of this to the N-H bond is covered previously
    if ((lcatom == "C" && hcatom == "N") ||
        (lcatom == "N" && hcatom == "C")) {
        return name + "CA_N_C_O";
    }

    // No special code for SG-SG
    return cresidue->name() + "_" + llatom + "_" + lcatom + "_" + hcatom + "_" +
           hhatom;
}

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    double bin_size = o.distance();

    std::unordered_map<std::thread::id, DihedralCounts> bins;

    auto worker = [bin_size, &bins](chemfiles::Frame complex,
                                    const std::string& /* unused */) {
        // Selection phase
        chemfiles::Frame protein_only;
        auto peptides =
            lemon::select_specific_residues(complex, lemon::common_peptides);

        if (peptides.size() == 0) {
            return;
        }

        lemon::separate_residues(complex, peptides, protein_only);
        protein_only.set_cell(complex.cell());

        const auto& dihedrals = protein_only.topology().dihedrals();

        for (const auto& dihedral : dihedrals) {
            std::string angle_name;
            auto improper_name = get_dihedral_name(protein_only, dihedral);

            auto theta = protein_only.dihedral(dihedral[0], dihedral[1],
                                               dihedral[2], dihedral[3]);

            int bin = static_cast<int>(std::floor(theta / bin_size));

            auto th = std::this_thread::get_id();
            BondDihedralBin sbin = {improper_name, bin};
            auto bin_iterator = bins[th].find(sbin);

            if (bin_iterator == bins[th].end()) {
                bins[th][sbin] = 1;
                continue;
            }

            ++(bin_iterator->second);
        }
    };

    auto p = o.work_dir();
    auto entries = o.entries();

    if (!boost::filesystem::is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    lemon::run_hadoop(worker, p);

    DihedralCounts sc_total;

    for (const auto& bin : bins) {
        for (const auto& sc : bin.second) {
            sc_total[sc.first] += sc.second;
        }
    }

    for (const auto& i : sc_total) {
        std::cout << i.first.first << "\t" << i.first.second * bin_size << "\t"
                  << i.second << "\n";
    }
}
