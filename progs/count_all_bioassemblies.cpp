#include <iostream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "benchmarker/entries.hpp"
#include "benchmarker/parse.hpp"
#include "benchmarker/run.hpp"
#include "benchmarker/select.hpp"

using namespace boost::filesystem;

typedef std::unordered_map<std::string, size_t> entry_to_asmbio;

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

    std::vector<entry_to_asmbio> resn_counts(ncpu);
    auto worker = [&resn_counts](const chemfiles::Frame& complex,
                                 const std::string& pdbid, size_t id) {

        auto result = benchmarker::count_bioassemblies(complex);
        resn_counts[id].emplace(pdbid, result);
    };

    current_path(p);
    benchmarker::call_multithreaded(worker, vec, ncpu, chun);

    for (const auto& iter : resn_counts) {
        for (const auto& iter2 : iter) {
            std::cout << iter2.first << "\t" << iter2.second << std::endl;
        }
    }
}
