#include <iostream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/count.hpp"
#include "lemon/entries.hpp"
#include "lemon/hadoop.hpp"
#include "lemon/options.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    std::unordered_map<std::thread::id, lemon::ResidueNameCount> resn_counts;
    auto worker = [&resn_counts](chemfiles::Frame complex,
                                 const std::string& /* unused */) {
        // Desired info is calculated directly, no pruning, output is done later
        auto th = std::this_thread::get_id();
        lemon::count_residues(complex, resn_counts[th]);
    };

    auto p = o.work_dir();
    auto threads = o.ncpu();

    if (!boost::filesystem::is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    try {
        lemon::run_hadoop(worker, p, threads);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }

    lemon::ResidueNameCount resn_total;

    for (const auto& resn_count : resn_counts) {
        resn_total += resn_count.second;
    }

    for (auto i : resn_total) {
        std::cout << i.first << "\t" << i.second << "\n";
    }
}
