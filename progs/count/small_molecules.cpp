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
        std::cout << lemon::count::print_residue_name_counts(pdbid, complex, smallm);
    };

    auto p = o.work_dir();
    auto threads = o.ncpu();

    try {
        lemon::run_parallel(worker, p, threads);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }
}
