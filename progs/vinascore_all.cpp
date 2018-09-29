#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/count.hpp"
#include "lemon/entries.hpp"
#include "lemon/hadoop.hpp"
#include "lemon/options.hpp"
#include "lemon/prune.hpp"
#include "lemon/score.hpp"
#include "lemon/select.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex,
                     const std::string& pdbid) {

        // Selection phase
        auto smallm = lemon::select_small_molecules(complex);

        // Pruning phase
        lemon::remove_identical_residues(complex, smallm);
        lemon::remove_cofactors(complex, smallm, lemon::common_cofactors);
        lemon::remove_cofactors(complex, smallm, lemon::linear_molecules);

        // Output phase
        const auto& residues = complex.topology().residues();
        std::set<size_t> proteins;
        for (size_t i = 0; i < complex.topology().residues().size(); ++i) {
            proteins.insert(i);
        }

        for (auto smallm_id : smallm) {
            auto prot_copy = proteins;

            lemon::keep_interactions(complex, smallm, prot_copy, 8.0);
            prot_copy.erase(smallm_id);

            auto vscore =
                lemon::xscore::vina_score(complex, smallm_id, prot_copy);

            std::stringstream ss;
            ss << pdbid << "\t" << residues[smallm_id].name() << "\t"
               << vscore.g1 << "\t" << vscore.g2 << "\t" << vscore.hydrogen
               << "\t" << vscore.hydrophobic << "\t" << vscore.rep << "\n";
            std::cout << ss.str();
        }
    };

    auto p = o.work_dir();
    auto entries = o.entries();

    if (!boost::filesystem::is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    lemon::run_hadoop(worker, p);
}
