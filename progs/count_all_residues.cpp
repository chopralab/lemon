#include <iostream>
#include <fstream>
#include <thread>
#include <iterator>

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

#include <chemfiles.hpp>

#include "benchmarker/parse.hpp"
#include "benchmarker/entries.hpp"

using namespace boost::filesystem;

int main(int argc, char *argv[]) {

    path entries(argc>1? argv[1] : "entries.idx");
    path p(argc>2? argv[2] : ".");
    size_t ncpu = argc>3? std::strtoul(argv[3], nullptr, 0) : 1;

    if(!is_regular_file(entries)) {
        std::cerr << "You must supply a valid entries file" << std::endl;
        return 1;
    }

    if(!is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    std::vector<std::array<char, 4>> vec;
    vec.reserve(142000); // increase as neccessisary

    std::ifstream input_file(entries.string());
    std::string header;
    std::string junk;
    std::getline(input_file, header);
    std::getline(input_file, junk);
    benchmarker::read_entry_file(input_file, vec);

    std::vector<benchmarker::ResidueNameCount> resn_counts(ncpu);

    typedef std::vector<std::array<char, 4>>::const_iterator iter;

    auto worker = [&p,&resn_counts] (iter begin, iter end, size_t id) {
        for(auto it = begin; it != end; ++it) {
        
            path entry = p / std::string(it->data(), 4);
            entry += std::string(".mmtf.gz");

            try {
                auto traj = chemfiles::Trajectory(entry.string());

                if (traj.nsteps() > 1) {
                    std::cerr << "File: " << entry << " contrains more than more model. Only considering the first" << std::endl;
                }

                auto complex = traj.read();
                benchmarker::retreive_residue_counts(complex, resn_counts[id]);
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
        }
    };

    std::vector<std::thread> threads(ncpu);
    const int grainsize = vec.size() / ncpu;

    auto work_iter = std::begin(vec);
    size_t id = 0;
    for(auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
        *it = std::thread(worker, work_iter, work_iter + grainsize, id++);
        work_iter += grainsize;
    }
    threads.back() = std::thread(worker, work_iter, std::end(vec), id++);

    for(auto&& i : threads) {
        i.join();
    }

    benchmarker::ResidueNameCount resn_total;

    for( const auto& resn_count : resn_counts) {
        resn_total += resn_count;
    }

    std::cout << resn_total;
}
