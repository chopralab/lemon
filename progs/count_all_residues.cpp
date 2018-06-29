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
    size_t complex_count = 0;

    for(auto& entry : boost::make_iterator_range(directory_iterator(p), {})) {
        try {
            auto traj = chemfiles::Trajectory(entry.path().string());

            if (traj.nsteps() > 1) {
                std::cerr << "File: " << entry << " contrains more than more model. Only considering the first" << std::endl;
            }

            auto complex = traj.read();
            benchmarker::retreive_residue_counts(complex, resn_count);
        } catch (const std::range_error& e) {
            std::cerr << "Odd residue name in " << entry << std::endl;
        } catch (const std::length_error& e) {
            std::cerr << "Long residue name in " << entry << std::endl;
        } catch (const chemfiles::FormatError& e) {
            std::cerr << "Unsupported file format for " << entry << ". Skipping" << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Unknown error: " << e.what() << " for " << entry << std::endl;
        } catch (...) {
            std::cerr << "Unknown error for " << entry << std::endl;
        }

        ++complex_count;

        if(complex_count % 1000 == 0) {
            std::cerr << "#";
        }
    }

    std::cerr << std::endl;
    std::cout << resn_count;
}
