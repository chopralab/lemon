#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

#include <chemfiles.hpp>
#include "benchmarker/parse.hpp"

using namespace boost::filesystem;

int main(int argc, char *argv[]) {
    path p(argc>1? argv[1] : ".");

    if(!is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 1;
    }

    benchmarker::ResidueNameCount resn_count;

    for(auto& entry : boost::make_iterator_range(directory_iterator(p), {})) {
        auto traj = chemfiles::Trajectory(entry.path().string());

        if (traj.nsteps() > 1) {
            std::cerr << "File: " << entry << " contrains more than more model. Only considering the first" << std::endl;
        }

        auto complex = traj.read();
        benchmarker::retreive_residue_counts(complex, resn_count);
    }

    std::cout << resn_count;
}
