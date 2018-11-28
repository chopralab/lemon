#include <iostream>
#include <sstream>

#include "lemon/count.hpp"
#include "lemon/entries.hpp"
#include "lemon/select.hpp"
#include "lemon/options.hpp"
#include "lemon/hadoop.hpp"
#include "lemon/parallel.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex,
                     const std::string& pdbid) {

        // Desired info is obtained directly
        auto result = lemon::count_bioassemblies(complex);

        // Custom output phase
        std::stringstream ss;
        ss << pdbid << " " << result << "\n";
        std::cout << ss.str();
    };

    auto p = o.work_dir();
    auto entries = o.entries();
    auto threads = o.ncpu();

    try {
        lemon::run_parallel(worker, p, threads);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }
}
