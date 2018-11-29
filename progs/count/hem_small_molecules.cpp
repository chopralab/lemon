#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    double distance = o.distance();

    auto worker = [distance](chemfiles::Frame complex,
                             const std::string& pdbid) {

        // Selection phase
        auto hemegs = lemon::select::specific_residues(
            complex, {"HEM", "HEA", "HEB", "HEC"});
        auto smallm = lemon::select::small_molecules(complex);

        // Pruning phase
        lemon::prune::identical_residues(complex, smallm);
        lemon::prune::cofactors(complex, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(complex, smallm, lemon::linear_molecules);

        lemon::prune::keep_interactions(complex, smallm, hemegs, distance);

        // Output phase
        lemon::count::print_residue_name_counts(std::cout, pdbid, complex, smallm);
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
