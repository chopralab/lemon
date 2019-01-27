#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex,
                     const std::string& pdbid) {

        // Selection phase
        std::list<size_t> smallm;
        if (lemon::select::small_molecules(complex, smallm) == 0) {
            return std::string("");
        }

        // Pruning phase
        lemon::prune::identical_residues(complex, smallm);
        lemon::prune::cofactors(complex, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(complex, smallm, lemon::common_fatty_acids);

        // Output phase
        return pdbid + lemon::count::print_residue_names(complex, smallm);
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
