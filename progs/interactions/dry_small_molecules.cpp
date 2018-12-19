#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto distance = 6.0;
    o.add_option("distance,d", distance,
				 "Largest distance between water and a small molecule.");
    o.parse_command_line(argc, argv);

    auto worker = [distance](chemfiles::Frame complex,
                             const std::string& pdbid) {

        // Selection phase
        auto waters = lemon::select::specific_residues(complex, {"HOH"});
        auto smallm = lemon::select::small_molecules(complex);

        // Pruning phase
        lemon::prune::identical_residues(complex, smallm);
        lemon::prune::cofactors(complex, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(complex, smallm, lemon::common_fatty_acids);

        lemon::prune::remove_interactions(complex, smallm, waters, distance);

        // Output phase
        return lemon::count::print_residue_name_counts(pdbid, complex, smallm);
    };

    return lemon::launch<lemon::print_combine>(o, worker, std::cout);
}
