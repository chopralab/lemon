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

        std::stringstream ss;
        ss << pdbid << "\t" << score << "\t" << rmsd << "\t" << aligned << "\n";
        std::cout << ss.str();
    };

    auto p = o.work_dir();
    auto entries = o.entries();
    auto threads = o.ncpu();

    try {
        lemon::run_parallel(worker, p, threads);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }
}
