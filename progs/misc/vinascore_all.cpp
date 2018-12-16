#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex,
                     const std::string& pdbid) {

        // Selection phase
        auto smallm = lemon::select::small_molecules(complex);

        // Pruning phase
        lemon::prune::identical_residues(complex, smallm);
        lemon::prune::cofactors(complex, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(complex, smallm, lemon::common_fatty_acids);

        // Output phase
        const auto& residues = complex.topology().residues();
        std::set<size_t> proteins;
        for (size_t i = 0; i < complex.topology().residues().size(); ++i) {
            proteins.insert(i);
        }

        for (auto smallm_id : smallm) {
            auto prot_copy = proteins;

            lemon::prune::keep_interactions(complex, smallm, prot_copy, 8.0);
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
    auto threads = o.ncpu();

    try {
        lemon::run_parallel(worker, p, threads);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }
}
