#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "benchmarker/entries.hpp"
#include "benchmarker/parse.hpp"
#include "benchmarker/run.hpp"
#include "benchmarker/select.hpp"

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

    auto worker = [](const chemfiles::Frame& complex,
                     const std::string& pdbid) {
        auto result = benchmarker::select_metal_ions(complex);

        if (result.empty()) {
            return;
        }

        const auto& residues = complex.topology().residues();

        std::unordered_map<std::string, size_t> metals;
        for (auto res_id : result) {
            auto atom_id = *(residues[res_id].begin());
            auto iter = metals.find(complex[atom_id].type());
            if (iter == metals.end()) {
                metals[complex[atom_id].type()] = 1;
                continue;
            }
            ++iter->second;
        }

        std::stringstream ss;
        ss << pdbid;
        for (const auto iter : metals) {
            ss << " " << iter.first << " " << iter.second;
        }
        ss << "\n";

        std::cout << ss.str();
    };

    current_path(p);
    benchmarker::call_multithreaded(worker, vec, ncpu, chun);
}
