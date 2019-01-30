#include <iostream>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame entry, const std::string&) {
        // Desired info is calculated directly, no pruning, output is done later
        lemon::ResidueNameCount rnc;
        lemon::count::residues(entry, rnc);
        return rnc;
    };

    lemon::ResidueNameCount resn_total;
    auto collector = lemon::map_combine<lemon::ResidueNameCount>(resn_total);
    lemon::launch(o, worker, collector);

    for (auto i : resn_total) {
        std::cout << i.first << "\t" << i.second << "\n";
    }
}

