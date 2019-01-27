#include <iostream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex, const std::string&) {

        lemon::ResidueNameCount rnc;

        // Selection phase
        chemfiles::Frame protein_only;
        auto peptides = lemon::select::peptides(complex);

        // Pruning phase
        if (peptides.size() == 0) {
            return rnc;
        }

        lemon::separate::residues(complex, peptides, protein_only);

        // Output phase
        lemon::count::residues(protein_only, rnc);
        return rnc;
    };

    lemon::ResidueNameCount resn_total;
    auto collector = lemon::map_combine<lemon::ResidueNameCount>(resn_total);
    lemon::launch(o, worker, collector);

    for (auto i : resn_total) {
        std::cout << i.first << "\t" << i.second << "\n";
    }
}
