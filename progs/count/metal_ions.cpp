#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame entry,
                     const std::string& pdbid) {

        // Selection phase
        std::list<size_t> metal_ids;
        lemon::select::metal_ions(entry, metal_ids);

        // No pruning, straight to out output phase
        return pdbid + lemon::count::print_residue_names(entry, metal_ids);
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
