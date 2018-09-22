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
typedef std::pair<std::string, int> BondImproperBin;
typedef std::map<BondImproperBin, size_t> ImproperCounts;

inline std::string get_improper_name(const chemfiles::Frame& complex,
                                     const chemfiles::Improper& improper) {
    const auto& atom1 = complex[improper[0]];
    const auto& atom2 = complex[improper[1]];
    const auto& atom3 = complex[improper[2]];
    const auto& atom4 = complex[improper[3]];

    const std::string& latom =
        std::min(atom1.name(), std::min(atom3.name(), atom4.name()));

    const std::string& matom =
        std::max(std::min(atom1.name(), atom3.name()),
                 std::min(std::max(atom1.name(), atom3.name()), atom4.name()));

    const std::string& hatom =
        std::max(atom1.name(), std::max(atom3.name(), atom4.name()));

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

    if (catom == "N") {
        return name + latom + "_" + catom + "_" + matom + "_" + hatom;
    }

    return cresidue->name() + "_" + latom + "_" + catom + "_" + matom + "_" +
           hatom;
}

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    double bin_size = o.distance();

    std::unordered_map<std::thread::id, ImproperCounts> bins;

    auto worker = [bin_size, &bins](chemfiles::Frame complex,
                                    const std::string& /* unused */) {
        // Selection phase
        chemfiles::Frame protein_only;
        auto peptides = lemon::select_peptides(complex);

        if (peptides.size() == 0) {
            return;
        }

        lemon::separate_residues(complex, peptides, protein_only);
        protein_only.set_cell(complex.cell());

        const auto& impropers = protein_only.topology().impropers();

        for (const auto& improper : impropers) {
            std::string angle_name;
            auto improper_name = get_improper_name(protein_only, improper);

            auto theta = protein_only.out_of_plane(improper[0], improper[1],
                                                   improper[2], improper[3]);

            int bin = static_cast<int>(std::floor(theta / bin_size));

            auto th = std::this_thread::get_id();
            BondImproperBin sbin = {improper_name, bin};
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

    ImproperCounts sc_total;

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
