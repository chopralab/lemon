#include <iostream>
#include <sstream>
#include <unordered_map>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"
#include "lemon/geometry.hpp"

// typedefs for binned data
typedef std::pair<std::string, int> BondDihedralBin;
typedef std::map<BondDihedralBin, size_t> DihedralCounts;

using lemon::geometry::protein::dihedral_name;

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto bin_size = 0.01;
    o.add_option("--bin_size,-b", bin_size, "Size of the dihedral bin.");
    o.parse_command_line(argc, argv);

    auto worker = [bin_size](chemfiles::Frame entry,
                             const std::string& pdbid) {
        DihedralCounts bins;
        
        // Selection phase
        chemfiles::Frame protein_only;
        std::list<size_t> peptides;

        if (lemon::select::specific_residues(entry, peptides,
                                             lemon::common_peptides) == 0) {
            return bins;
        }

        lemon::separate::residues(entry, peptides, protein_only);
        protein_only.set_cell(entry.cell());

        const auto& dihedrals = protein_only.topology().dihedrals();

        for (const auto& dihedral : dihedrals) {
            std::string dihedralnm;
            try {
                dihedralnm = dihedral_name(protein_only, dihedral);
            } catch (const lemon::geometry::geometry_error& e) {
                auto msg = pdbid + ": " + e.what() + '\n';
                std::cerr << msg;
            }

            auto theta = protein_only.dihedral(dihedral[0], dihedral[1],
                                               dihedral[2], dihedral[3]);

            int bin = static_cast<int>(std::floor(theta / bin_size));

            BondDihedralBin sbin = {dihedralnm, bin};
            auto bin_iterator = bins.find(sbin);

            if (bin_iterator == bins.end()) {
                bins[sbin] = 1;
                continue;
            }

            ++(bin_iterator->second);
        }

        return bins;
    };

    DihedralCounts sc_total;
    auto collector = lemon::map_combine<DihedralCounts>(sc_total);
    lemon::launch(o, worker, collector);

    for (const auto& i : sc_total) {
        std::cout << i.first.first << "\t"
                  << static_cast<double>(i.first.second) * bin_size << "\t"
                  << i.second << "\n";
    }

    return 0;
}
