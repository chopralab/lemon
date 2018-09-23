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

    // Let's keep the aminoacid group all the same (except proline)
    if (latom == "CA" || hatom == "CA") {
        return name + latom + "_" + hatom;
    }

    return residue1->name() + "_" + latom + "_" + hatom;
}

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    double bin_size = o.distance();

    std::unordered_map<std::thread::id, StretchCounts> bins;

    auto worker = [bin_size, &bins](chemfiles::Frame complex,
                                    const std::string& /* unused */) {
        // Selection phase
        chemfiles::Frame protein_only;
        auto peptides = lemon::select_peptides(complex);

        if (peptides.size() == 0) {
            return;
        }

        lemon::separate_residues(complex, peptides, protein_only);

        const auto& bonds = protein_only.topology().bonds();

        for (const auto& bond : bonds) {
            std::string bond_name = get_bond_name(protein_only, bond);

            auto distance = protein_only.distance(bond[0], bond[1]);
            size_t bin = static_cast<size_t>(std::floor(distance / bin_size));

            auto th = std::this_thread::get_id();
            BondStretchBin sbin = {bond_name, bin};
            auto bin_iterator = bins[th].find(sbin);

            if (bin_iterator == bins[th].end()) {
                bins[th][sbin] = 1;
                continue;
            }

            ++(bin_iterator->second);
        }
    };

    auto p = o.work_dir();
    auto ncpu = o.npu();
    auto entries = o.entries();

    if (!boost::filesystem::is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    lemon::run_hadoop(worker, p, ncpu);

    StretchCounts sc_total;

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
