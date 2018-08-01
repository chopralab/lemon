#include <iostream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/entries.hpp"
#include "lemon/count.hpp"
#include "lemon/archive_run.hpp"
#include "lemon/options.hpp"
#include "lemon/hadoop.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    std::unordered_map<std::thread::id, lemon::ResidueNameCount>
        resn_counts;
    auto worker = [&resn_counts](const chemfiles::Frame& complex,
                                 const std::string& /* unused */) {

        // Desired info is calculated directly, no pruning, output is done
        auto th = std::this_thread::get_id();
        lemon::count_residues(complex, resn_counts[th]);
    };

    auto p = o.work_dir();
    auto ncpu = o.npu();
    auto entries = o.entries();

    if (!boost::filesystem::is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    if (!entries.empty()) {

        if (!boost::filesystem::is_regular_file(entries)) {
            std::cerr << "You must supply a valid entries file" << std::endl;
            return 1;
        }

        lemon::PDBIDVec vec;
        lemon::read_entry_file(entries, vec);
        boost::filesystem::current_path(p);
        lemon::run_archive(worker, vec, ncpu, 1);
    } else {
        lemon::run_hadoop(worker, p, ncpu);
    }

    lemon::ResidueNameCount resn_total;

    for (const auto& resn_count : resn_counts) {
        resn_total += resn_count.second;
    }

    for (auto i : resn_total) {
        std::cout << i.first << "\t" << i.second << "\n";
    }
}
