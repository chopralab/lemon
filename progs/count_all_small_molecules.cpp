#include <iostream>
#include <fstream>
#include <thread>
#include <iterator>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "benchmarker/parse.hpp"
#include "benchmarker/entries.hpp"
#include "benchmarker/select.hpp"
#include "benchmarker/prune.hpp"

using namespace boost::filesystem;

typedef std::unordered_map<std::string, std::unordered_map<std::string, size_t>> entry_to_small_molecule;

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

    std::vector<entry_to_small_molecule> resn_counts(ncpu);

    typedef std::vector<std::array<char, 4>>::const_iterator iter;

    auto worker = [&p,&resn_counts] (iter begin, iter end, size_t id) {
        for(auto it = begin; it != end; ++it) {
            path entry = p / std::string(it->data(), 1) / std::string(it->data() + 1, 1) / std::string(it->data(), 4);
            entry += std::string(".mmtf.gz");

            try {
                auto traj = chemfiles::Trajectory(entry.string());

                if (traj.nsteps() > 1) {
                    std::cerr << "File: " << entry << " contrains more than more model. Only considering the first" << std::endl;
                }

                auto complex = traj.read();
                auto result = benchmarker::select_small_molecule(complex);
                benchmarker::remove_identical_residues(complex, result);
                benchmarker::remove_common_cofactors(complex, result);
                std::unordered_map<std::string, size_t> small_molecules;

                const auto& residues = complex.topology().residues();

                for (auto atom_id : result) {
                    auto iter = small_molecules.find(residues[atom_id].name());
                    if (iter == small_molecules.end()) {
                        small_molecules[residues[atom_id].name()] = 1;
                        continue;
                    }
                    ++iter->second;
                }
                resn_counts[id].emplace(std::string(it->data(), 4), small_molecules);
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

    for (const auto& iter : resn_counts) {
        for (const auto& iter2 : iter) {

            if (iter2.second.size() == 0) {
                continue;
            }

            std::cout << iter2.first;

            for (const auto iter3 : iter2.second) {
                std::cout << " " << iter3.first << " " << iter3.second;
            }

            std::cout << std::endl;
        }
    }
}
