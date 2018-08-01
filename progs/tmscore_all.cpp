#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/count.hpp"
#include "lemon/entries.hpp"
#include "lemon/prune.hpp"
#include "lemon/archive_run.hpp"
#include "lemon/select.hpp"
#include "lemon/structure.hpp"

using namespace boost::filesystem;

int main(int argc, char* argv[]) {
    path entries(argc > 1 ? argv[1] : "entries.idx");
    path p(argc > 2 ? argv[2] : ".");
    path reference(argc > 3 ? argv[3] : "reference.mmtf.gz");
    size_t ncpu = argc > 4 ? std::strtoul(argv[4], nullptr, 0) : 1;
    size_t chun = argc > 5 ? std::strtoul(argv[5], nullptr, 0) : 1;

    if (!is_regular_file(entries)) {
        std::cerr << "You must supply a valid entries file" << std::endl;
        return 1;
    }

    if (!is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    if (!is_regular_file(reference)) {
        std::cerr << "You must supply a valid reference file" << std::endl;
        return 3;
    }

    chemfiles::Trajectory traj(reference.string());
    chemfiles::Frame native = traj.read();

    lemon::PDBIDVec vec;
    lemon::read_entry_file(entries.string(), vec);

    auto worker = [&native](const chemfiles::Frame& complex,
                            const std::string& pdbid) {

        std::vector<chemfiles::Vector3D> junk;

        double score, rmsd;
        size_t aligned;
        std::tie(score, rmsd, aligned) = lemon::TMscore(complex, native, junk);

        std::stringstream ss;
        ss << pdbid << "\t" << score << "\t" << rmsd << "\t" << aligned << "\n";
        std::cout << ss.str();
    };

    current_path(p);
    lemon::run_archive(worker, vec, ncpu, chun);
}
