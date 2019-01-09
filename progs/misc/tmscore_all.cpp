#include <iostream>
#include <sstream>

#include "lemon/lemon.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto reference = std::string("reference.pdb");
    o.add_option("reference,r", reference, "Protein or DNA to align to.");
    o.parse_command_line(argc, argv);

    chemfiles::Trajectory traj(reference);
    chemfiles::Frame native = traj.read();

    auto worker = [&native](chemfiles::Frame complex,
                            const std::string& pdbid) {

        std::vector<chemfiles::Vector3D> junk;

        double score, rmsd;
        size_t aligned;
        std::tie(score, rmsd, aligned) = lemon::tmalign::TMscore(complex, native, junk);

        return pdbid + "\t" + std::to_string(score) + "\t" + std::to_string(rmsd) + "\t" + std::to_string(aligned) + "\n";
    };

    return lemon::launch<lemon::print_combine>(o, worker, std::cout);
}
