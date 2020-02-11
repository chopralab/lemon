#include <iostream>
#include <sstream>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto distance = lemon::prune::DEFAULT_DISTANCE;
    o.add_option("--distance,-d", distance,
                 "Largest distance between a SAM group and a small molecule.");
    o.parse_command_line(argc, argv);

    auto worker = [distance](const chemfiles::Frame& entry,
                             const std::string& pdbid) {

        // Selection phase
        auto sam = lemon::select::specific_residues(entry, {"SAM"});
        auto smallm = lemon::select::small_molecules(entry);

        // Pruning phase
        lemon::prune::identical_residues(entry, smallm);
        lemon::prune::cofactors(entry, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(entry, smallm, lemon::common_fatty_acids);

        lemon::prune::keep_interactions(entry, smallm, sam, distance);

        // Output phase
        return pdbid + lemon::count::print_residue_names(entry, smallm);
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
