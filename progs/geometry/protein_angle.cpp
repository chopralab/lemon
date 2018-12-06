#include <iostream>
#include <sstream>
#include <unordered_map>

#include "lemon/lemon.hpp"

// typedefs for binned data
typedef std::pair<std::string, size_t> BondAngleBin;
typedef std::map<BondAngleBin, size_t> AngleCounts;

inline std::string get_angle_name(const chemfiles::Frame& complex,
                                  const chemfiles::Angle& angle) {
    const auto& atom1 = complex[angle[0]];
    const auto& atom2 = complex[angle[1]];
    const auto& atom3 = complex[angle[2]];

    if (atom1.type() == "H") {
        return atom3.type() + "_" + atom2.type() + "_H";
    }

    if (atom3.type() == "H") {
        return atom1.type() + "_" + atom2.type() + "_H";
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
                return "N_C_O";
            }
            if (latom == "OXT" || hatom == "OXT") {
                return "O_C_OXT";
            }
            // This should not be possible....
        }

        // Maybe proline logic should be added here?
        return "CA_C_N";
    }

    if ((latom == "C" || hatom == "C") && catom == "CA") {
        return "C_CA_N";
    }

    if ((latom == "C" || hatom == "C") && catom == "N") {
        if (latom == "CA" || hatom == "CA") {
            return "C_N_CA";
        }
        // Proline!
        if (latom == "CD" || hatom == "CD") {
            return "PRO_C_N_CD";
        }
        return "C_N_H";  // Hydrogens???
    }

    // Check for sulfide bridges
    if (catom == "SG" && (latom == "SG" || hatom == "SG")) {
        return "CB_SG_SG";
    }

    const auto& residue1 = complex.topology().residue_for_atom(angle[0]);
    const auto& residue3 = complex.topology().residue_for_atom(angle[2]);

    if (residue1 != residue3) {
        std::cerr << "Unhandled inter-residue angle: " << residue1->name()
                  << " " << atom1.name() << " " << residue3->name() << " "
                  << atom3.name() << "\n";
        return "err";
    }

    return residue1->name() + "_" + latom + "_" + catom + "_" + hatom;
}

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    double bin_size = o.distance();

    auto worker = [bin_size](chemfiles::Frame complex, const std::string&) {
        AngleCounts bins;

        // Selection phase
        chemfiles::Frame protein_only;
        auto peptides =
            lemon::select::specific_residues(complex, lemon::common_peptides);

        if (peptides.size() == 0) {
            return bins;
        }

        lemon::separate::residues(complex, peptides, protein_only);

        const auto& angles = protein_only.topology().angles();

        for (const auto& angle : angles) {
            std::string angle_name;
            angle_name = get_angle_name(protein_only, angle);

            auto theta = protein_only.angle(angle[0], angle[1], angle[2]);
            size_t bin = static_cast<size_t>(std::floor(theta / bin_size));

            BondAngleBin sbin = {angle_name, bin};
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

    lemon::map_combine<AngleCounts> combiner;
    AngleCounts sc_total;

    try {
        lemon::run_parallel(worker, combiner, p, sc_total, threads);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }

    for (const auto& i : sc_total) {
        std::cout << i.first.first << "\t"
                  << static_cast<double>(i.first.second) * bin_size << "\t"
                  << i.second << "\n";
    }
}
