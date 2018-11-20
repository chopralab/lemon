#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/entries.hpp"
#include "lemon/count.hpp"
#include "lemon/prune.hpp"
#include "lemon/select.hpp"
#include "lemon/options.hpp"
#include "lemon/hadoop.hpp"

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
        lemon::print_residue_name_counts(std::cout, pdbid, complex, smallm);
    };

    auto p = o.work_dir();
    auto threads = o.ncpu();

    try {
        lemon::run_hadoop(worker, p, threads);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }
}
