#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex,
                     const std::string& pdbid) {

        // Selection phase
        auto result = lemon::select::metal_ions(complex);

        // No pruning, straight to out output phase
        return lemon::count::print_residue_name_counts(pdbid, complex, result);
    };

    return lemon::launch<lemon::print_combine>(o, worker, std::cout);
}
