#include <iostream>
#include <sstream>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](const chemfiles::Frame& entry,
                     const std::string& pdbid) -> std::string {

        // Selection phase
        auto metal_ids = lemon::select::metal_ions(entry);

        // No pruning, straight to out output phase
        return pdbid + lemon::count::print_residue_names(entry, metal_ids);
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
