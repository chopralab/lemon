#include <iostream>
#include <sstream>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](const chemfiles::Frame& entry,
                     const std::string& pdbid) -> std::string {

        // Selection phase
        auto smallm = lemon::select::small_molecules(entry);
        if (smallm.empty()) {
            return std::string("");
        }

        // Pruning phase
        lemon::prune::identical_residues(entry, smallm);
        lemon::prune::cofactors(entry, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(entry, smallm, lemon::common_fatty_acids);

        // Output phase
        return pdbid + lemon::count::print_residue_names(entry, smallm);
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
