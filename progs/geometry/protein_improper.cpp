#include <iostream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "lemon/lemon.hpp"

// typedefs for binned data
typedef std::pair<std::string, int> BondImproperBin;
typedef std::map<BondImproperBin, size_t> ImproperCounts;

using lemon::geometry::protein::improper_name;

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto bin_size = 0.01;
    o.add_option("bin_size,b", bin_size, "Size of the improper-dihedral bin.");
    o.parse_command_line(argc, argv);

    auto worker = [bin_size](chemfiles::Frame complex,
                             const std::string& pdbid) {
        ImproperCounts bins;

        // Selection phase
        chemfiles::Frame protein_only;
        std::list<size_t> peptides;

        if (lemon::select::specific_residues(complex, peptides,
                                             lemon::common_peptides) == 0) {
            return bins;
        }

        lemon::separate::residues(complex, peptides, protein_only);
        protein_only.set_cell(complex.cell());

        const auto& impropers = protein_only.topology().impropers();

        for (const auto& improper : impropers) {
            std::string impropernm;
            try {
                impropernm = improper_name(protein_only, improper);
            } catch (const lemon::geometry::geometry_error& e) {
                auto msg = pdbid + ": " + e.what() + '\n';
                std::cerr << msg;
            }

            auto theta = protein_only.out_of_plane(improper[0], improper[1],
                                                   improper[2], improper[3]);

            int bin = static_cast<int>(std::floor(theta / bin_size));

            BondImproperBin sbin = {impropernm, bin};
            auto bin_iterator = bins.find(sbin);

            if (bin_iterator == bins.end()) {
                bins[sbin] = 1;
                continue;
            }

            ++(bin_iterator->second);
        }

        return bins;
    };

    ImproperCounts sc_total;
    lemon::launch<lemon::map_combine>(o, worker, sc_total);

    for (const auto& i : sc_total) {
        std::cout << i.first.first << "\t"
                  << static_cast<double>(i.first.second) * bin_size << "\t"
                  << i.second << "\n";
    }

    return 0;
}
