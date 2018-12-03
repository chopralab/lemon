#include <iostream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "lemon/lemon.hpp"

// typedefs for binned data
typedef std::pair<std::string, int> BondImproperBin;
typedef std::map<BondImproperBin, size_t> ImproperCounts;

inline const std::string& min_str(const std::string& s1, const std::string& s2) {
    return (s1 < s2)? s1 : s2;
}

inline const std::string& max_str(const std::string& s1, const std::string& s2) {
    return (s1 > s2)? s1 : s2;
}

inline std::string get_improper_name(const chemfiles::Frame& complex,
                                     const chemfiles::Improper& improper) {
    const auto& atom1 = complex[improper[0]];
    const auto& atom2 = complex[improper[1]];
    const auto& atom3 = complex[improper[2]];
    const auto& atom4 = complex[improper[3]];

    const std::string& latom =
        min_str(atom1.name(), min_str(atom3.name(), atom4.name()));

    const std::string& matom =
        max_str(min_str(atom1.name(), atom3.name()),
                min_str(max_str(atom1.name(), atom3.name()), atom4.name()));

    const std::string& hatom =
        max_str(atom1.name(), max_str(atom3.name(), atom4.name()));

    const std::string& catom = atom2.name();

    if (atom1.type() == "H") {
        return atom3.type() + "_" + atom2.type() + "_" + atom4.type() + "_H";
    }

    if (atom3.type() == "H") {
        return atom1.type() + "_" + atom2.type() + "_" + atom4.type() + "_H";
    }

    if (atom4.type() == "H") {
        return atom1.type() + "_" + atom2.type() + "_" + atom3.type() + "_H";
    }

    // The carbonyl angles should all be the same
    if (catom == "C") {
        // It must be the peptide bond!
        // as CA_C_N_H is already handled
        return "CA_C_N_O";
    }

    const auto& cresidue = complex.topology().residue_for_atom(improper[1]);

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

    if (catom == "N") {
        return name + latom + "_" + catom + "_" + matom + "_" + hatom;
    }

    return cresidue->name() + "_" + latom + "_" + catom + "_" + matom + "_" +
           hatom;
}

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    double bin_size = o.distance();

    auto worker = [bin_size](chemfiles::Frame complex, const std::string&) {
        ImproperCounts bins;

        // Selection phase
        chemfiles::Frame protein_only;
        auto peptides =
            lemon::select::specific_residues(complex, lemon::common_peptides);

        if (peptides.size() == 0) {
            return bins;
        }

        lemon::separate::residues(complex, peptides, protein_only);
        protein_only.set_cell(complex.cell());

        const auto& impropers = protein_only.topology().impropers();

        for (const auto& improper : impropers) {
            std::string angle_name;
            auto improper_name = get_improper_name(protein_only, improper);

            auto theta = protein_only.out_of_plane(improper[0], improper[1],
                                                   improper[2], improper[3]);

            int bin = static_cast<int>(std::floor(theta / bin_size));

            BondImproperBin sbin = {improper_name, bin};
            auto bin_iterator = bins.find(sbin);

            if (bin_iterator == bins.end()) {
                bins[sbin] = 1;
                continue;
            }

            ++(bin_iterator->second);
        }

        return bins;
    };

    auto p = o.work_dir();
    auto entries = o.entries();
    auto threads = o.ncpu();

    lemon::map_combine<ImproperCounts> combiner;
    ImproperCounts sc_total;

    try {
        lemon::run_parallel(worker, combiner, p, sc_total, threads);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }

    for (const auto& i : sc_total) {
        std::cout << i.first.first << "\t" << i.first.second * bin_size << "\t"
                  << i.second << "\n";
    }
}
