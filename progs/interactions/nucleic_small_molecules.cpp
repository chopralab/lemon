#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto distance = 6.0;
    o.add_option("distance,d", distance,
                 "Largest distance between a nucleic-acid and a small molecule.");
    o.parse_command_line(argc, argv);

    auto worker = [distance](chemfiles::Frame complex,
                             const std::string& pdbid) {

        // Selection phase
        auto nucleic_acids = lemon::select::nucleic_acids(complex);
        auto smallm = lemon::select::small_molecules(complex);

        // Pruning phase
        lemon::prune::identical_residues(complex, smallm);
        lemon::prune::cofactors(complex, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(complex, smallm, lemon::common_fatty_acids);
        lemon::prune::keep_interactions(complex, smallm, nucleic_acids, distance);

        // Output phase
        return lemon::count::print_residue_name_counts(pdbid, complex, smallm);
    };

    return lemon::launch<lemon::print_combine>(o, worker, std::cout);
}
