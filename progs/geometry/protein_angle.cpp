#include <iostream>
#include <sstream>
#include <unordered_map>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"
#include "lemon/geometry.hpp"

// typedefs for binned data
using BondAngleBin = std::pair<std::string, size_t>;
using AngleCounts  = std::map<BondAngleBin, size_t>;

using lemon::geometry::protein::angle_name;

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto bin_size = 0.01;
    o.add_option("--bin_size,-b", bin_size, "Size of the angle bin.");
    o.parse_command_line(argc, argv);

    auto worker = [bin_size](chemfiles::Frame entry,
                             const std::string& pdbid) {
        AngleCounts bins;

        // Selection phase
        chemfiles::Frame protein_only;
        std::list<size_t> peptides;

        if (lemon::select::specific_residues(entry, peptides,
                                             lemon::common_peptides) == 0) {
            return bins;
        }

        lemon::separate::residues(entry, peptides, protein_only);

        const auto& angles = protein_only.topology().angles();

        for (const auto& angle : angles) {
            std::string anglenm;
            try {
                anglenm = angle_name(protein_only, angle);
            } catch (const lemon::geometry::geometry_error& e) {
                auto msg = pdbid + ": " + e.what() + '\n';
                std::cerr << msg;
            }

            auto theta = protein_only.angle(angle[0], angle[1], angle[2]);
            auto bin = static_cast<size_t>(std::floor(theta / bin_size));

            BondAngleBin sbin = {anglenm, bin};
            auto bin_iterator = bins.find(sbin);

            if (bin_iterator == bins.end()) {
                bins[sbin] = 1;
                continue;
            }

            ++(bin_iterator->second);
        }

        return bins;
    };

    AngleCounts sc_total;
    auto collector = lemon::map_combine<AngleCounts>(sc_total);
    lemon::launch(o, worker, collector);

    for (const auto& i : sc_total) {
        std::cout << i.first.first << "\t"
                  << static_cast<double>(i.first.second) * bin_size << "\t"
                  << i.second << "\n";
    }

    return 0;
}
