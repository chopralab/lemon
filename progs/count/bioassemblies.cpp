#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex, const std::string& pdbid) {

        // Desired info is obtained directly
        auto result = lemon::count::bioassemblies(complex);

        // Output phase
        return pdbid + " " + std::to_string(result) + "\n";
    };

    return lemon::launch<lemon::print_combine>(o, worker, std::cout);
}
