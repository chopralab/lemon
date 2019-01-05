#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex,
                     const std::string& pdbid) {

        // Selection phase
        std::list<size_t> metal_ids;
        lemon::select::metal_ions(complex, metal_ids);

        // No pruning, straight to out output phase
        return lemon::count::print_residue_name_counts(pdbid, complex, metal_ids);
    };

    return lemon::launch<lemon::print_combine>(o, worker, std::cout);
}
