#include <iostream>
#include <sstream>
#include <unordered_map>

#include "lemon/lemon.hpp"

// typedefs for binned data
typedef std::pair<std::string, size_t> BondStretchBin;
typedef std::map<BondStretchBin, size_t> StretchCounts;

inline std::string get_bond_name(const chemfiles::Frame& complex,
                                 const chemfiles::Bond& bond) {
    const auto& atom1 = complex[bond[0]];
    const auto& atom2 = complex[bond[1]];

    if (atom1.type() == "H") {
        return atom2.type() + "_H";
    }

    if (atom2.type() == "H") {
        return atom1.type() + "_H";
    }

    const std::string latom =
        atom1.name() < atom2.name() ? atom1.name() : atom2.name();

    const std::string hatom =
        atom1.name() > atom2.name() ? atom1.name() : atom2.name();

    // The carbonyl bond should all be the same
    if (latom == "C" && hatom == "O") {
        return "C_O";
    }

    const auto& residue1 = complex.topology().residue_for_atom(bond[0]);
    const auto& residue2 = complex.topology().residue_for_atom(bond[1]);

    if (residue1 != residue2) {
        if (latom == "C" && hatom == "N") {
            return "peptide_bond";
        }

        if (latom == "SG" && hatom == "SG") {
            return "SG_SG";
        }

        std::cerr << "Unhandled inter-residue bond: " << residue1->name() << " "
                  << atom1.name() << " " << residue2->name() << " "
                  << atom2.name() << "\n";
        return "err";
    }

    std::string name;
    // Proline is ... special
    if (residue1->name() == "PRO") {
        name = "PRO_";
    }

    // Same goes for hydroxyproline
    if (residue1->name() == "HYP") {
        name = "HYP_";
    }

    // and pyroglutamic acid
    if (residue1->name() == "PCA") {
        name = "PCA_";
    }

    // Let's keep the aminoacid group all the same (except proline)
    if (latom == "CA" || hatom == "CA") {
        return name + latom + "_" + hatom;
    }

    return residue1->name() + "_" + latom + "_" + hatom;
}

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto bin_size = 0.01;
    o.add_option("bin_size,b", bin_size, "Size of the length(stretch) bin.");
    o.parse_command_line(argc, argv);

    auto worker = [bin_size](chemfiles::Frame complex, const std::string&) {
        StretchCounts bins;

        // Selection phase
        chemfiles::Frame protein_only;
        std::list<size_t> peptides;

        if (lemon::select::specific_residues(complex, peptides,
                                             lemon::common_peptides) == 0) {
            return bins;
        }

        lemon::separate::residues(complex, peptides, protein_only);

        const auto& bonds = protein_only.topology().bonds();

        for (const auto& bond : bonds) {
            std::string bond_name = get_bond_name(protein_only, bond);

            auto distance = protein_only.distance(bond[0], bond[1]);
            size_t bin = static_cast<size_t>(std::floor(distance / bin_size));

            BondStretchBin sbin = {bond_name, bin};
            auto bin_iterator = bins.find(sbin);

            if (bin_iterator == bins.end()) {
                bins[sbin] = 1;
                continue;
            }

            ++(bin_iterator->second);
        }

        return bins;
    };

    StretchCounts sc_total;
    lemon::launch<lemon::map_combine>(o, worker, sc_total);


    for (const auto& i : sc_total) {
        std::cout << i.first.first << "\t"
                  << static_cast<double>(i.first.second) * bin_size << "\t"
                  << i.second << "\n";
    }

    return 0;
}
