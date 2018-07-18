#include <iostream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "benchmarker/entries.hpp"
#include "benchmarker/parse.hpp"
#include "benchmarker/run.hpp"

using namespace boost::filesystem;

int main(int argc, char* argv[]) {
    path entries(argc > 1 ? argv[1] : "entries.idx");
    path p(argc > 2 ? argv[2] : ".");
    size_t ncpu = argc > 3 ? std::strtoul(argv[3], nullptr, 0) : 1;
    size_t chun = argc > 4 ? std::strtoul(argv[4], nullptr, 0) : 1;

    if (!is_regular_file(entries)) {
        std::cerr << "You must supply a valid entries file" << std::endl;
        return 1;
    }

    if (!is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    std::vector<std::array<char, 4>> vec;
    benchmarker::read_entry_file(entries.string(), vec);

    std::vector<benchmarker::ResidueNameCount> resn_counts(ncpu);
    auto worker = [&resn_counts](const chemfiles::Frame& complex,
                                 const std::string& pdbid, size_t id) {
        benchmarker::retreive_residue_counts(complex, resn_counts[id]);
    };

    current_path(p);
    benchmarker::call_multithreaded(worker, vec, ncpu, chun);

    benchmarker::ResidueNameCount resn_total;

    for (const auto& resn_count : resn_counts) {
        resn_total += resn_count;
    }

    std::cout << resn_total;
}
