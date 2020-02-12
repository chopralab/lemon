#include <iostream>
#include <sstream>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](const chemfiles::Frame& entry,
                     const std::string& pdbid) -> std::string {

        // Desired info is obtained directly 
        auto result = lemon::count::atom_property(entry, "altloc");

        // Output phase
        return pdbid + " " + std::to_string(result) + "\n";
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
