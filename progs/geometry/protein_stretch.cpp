#include <iostream>
#include <sstream>
#include <unordered_map>

#include "lemon/lemon.hpp"
#include "lemon/geometry.hpp"

// typedefs for binned data
typedef std::pair<std::string, size_t> BondStretchBin;
typedef std::map<BondStretchBin, size_t> StretchCounts;

using lemon::geometry::protein::bond_name;

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto bin_size = 0.01;
    o.add_option("bin_size,b", bin_size, "Size of the length(stretch) bin.");
    o.parse_command_line(argc, argv);

    auto worker = [bin_size](chemfiles::Frame complex,
                             const std::string& pdbid) {
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
            std::string bondnm;
            try {
                bondnm = bond_name(protein_only, bond);
            } catch (const lemon::geometry::geometry_error& e) {
                auto msg = pdbid + ": " + e.what() + '\n';
                std::cerr << msg;
            }
            auto distance = protein_only.distance(bond[0], bond[1]);
            size_t bin = static_cast<size_t>(std::floor(distance / bin_size));

            BondStretchBin sbin = {bondnm, bin};
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
    auto collector = lemon::map_combine<StretchCounts>(sc_total);
    lemon::launch(o, worker, collector);

    for (const auto& i : sc_total) {
        std::cout << i.first.first << "\t"
                  << static_cast<double>(i.first.second) * bin_size << "\t"
                  << i.second << "\n";
    }

    return 0;
}
