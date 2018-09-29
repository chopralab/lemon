#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/entries.hpp"
#include "lemon/count.hpp"
#include "lemon/select.hpp"
#include "lemon/options.hpp"
#include "lemon/hadoop.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    auto worker = [](chemfiles::Frame complex, const std::string& pdbid) {

        // Desired info is obtained directly 
        auto result = lemon::count_altloc(complex);

        // Custom output phase
        std::stringstream ss;
        ss << pdbid << " " << result << "\n";
        std::cout << ss.str();
    };

    auto p = o.work_dir();

    if (!boost::filesystem::is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    lemon::run_hadoop(worker, p);
}
