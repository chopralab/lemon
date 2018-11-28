#include <iostream>
#include <sstream>

#include "lemon/count.hpp"
#include "lemon/select.hpp"
#include "lemon/options.hpp"
#include "lemon/hadoop.hpp"
#include "lemon/parallel.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex,
                     const std::string& pdbid) {

        // Selection phase
        auto result = lemon::select_metal_ions(complex);

        // No pruning, straight to out output phase
        lemon::print_residue_name_counts(std::cout, pdbid, complex, result);
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
