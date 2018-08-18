#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/archive_run.hpp"
#include "lemon/count.hpp"
#include "lemon/entries.hpp"
#include "lemon/hadoop.hpp"
#include "lemon/options.hpp"
#include "lemon/prune.hpp"
#include "lemon/select.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o(argc, argv);

    double distance = o.distance();

    auto worker = [distance](chemfiles::Frame complex,
                             const std::string& pdbid) {

        // Selection phase
        auto sam = lemon::select_specific_residues(complex, {"SAM"});
        auto smallm = lemon::select_small_molecules(complex);

        // Pruning phase
        lemon::remove_identical_residues(complex, smallm);
        lemon::remove_cofactors(complex, smallm, lemon::common_cofactors);
        lemon::remove_cofactors(complex, smallm, lemon::linear_molecules);

        lemon::keep_interactions(complex, smallm, sam, distance);

        // Output phase
        lemon::print_residue_name_counts(std::cout, pdbid, complex, smallm);
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
}
