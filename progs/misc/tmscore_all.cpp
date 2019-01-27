#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"
#include "lemon/tmalign.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto reference = std::string("reference.pdb");
    o.add_option("reference,r", reference, "Protein or DNA to align to.");
    o.parse_command_line(argc, argv);

    chemfiles::Trajectory traj(reference);
    chemfiles::Frame native = traj.read();

    auto worker = [&native](chemfiles::Frame entry,
                            const std::string& pdbid) {

        std::vector<chemfiles::Vector3D> junk;

        auto tm = lemon::tmalign::TMscore(entry, native, junk);

        return pdbid + "\t" +
               std::to_string(tm.score) + "\t" +
               std::to_string(tm.rmsd) + "\t" +
               std::to_string(tm.aligned) + "\n";
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
