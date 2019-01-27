#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex, const std::string& pdbid) {

        // Desired info is obtained directly
        auto result = lemon::count::residue_property(complex, "assembly");

        // Output phase
        return pdbid + " " + std::to_string(result) + "\n";
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
