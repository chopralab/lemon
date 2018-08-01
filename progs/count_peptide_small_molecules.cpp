#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/entries.hpp"
#include "lemon/count.hpp"
#include "lemon/prune.hpp"
#include "lemon/archive_run.hpp"
#include "lemon/select.hpp"

using namespace boost::filesystem;

typedef std::unordered_map<std::string, std::unordered_map<std::string, size_t>>
    entry_to_small_molecule;

int main(int argc, char* argv[]) {
    path entries(argc > 1 ? argv[1] : "entries.idx");
    path p(argc > 2 ? argv[2] : ".");
    double dist_cutoff = argc > 3 ? std::strtod(argv[3], nullptr) : 6.0;
    size_t ncpu = argc > 4 ? std::strtoul(argv[4], nullptr, 0) : 1;
    size_t chun = argc > 5 ? std::strtoul(argv[5], nullptr, 0) : 1;

    if (!is_regular_file(entries)) {
        std::cerr << "You must supply a valid entries file" << std::endl;
        return 1;
    }

    if (!is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    lemon::PDBIDVec vec;
    lemon::read_entry_file(entries.string(), vec);

    auto worker = [dist_cutoff](const chemfiles::Frame& complex,
                                const std::string& pdbid) {

        // Selection phase
        auto peptides = lemon::select_peptides(complex);
        auto smallm = lemon::select_small_molecules(complex);

        // Pruning phase
        lemon::remove_identical_residues(complex, smallm);
        lemon::remove_cofactors(complex, smallm, lemon::common_cofactors);
        lemon::remove_cofactors(complex, smallm, lemon::linear_molecules);
        lemon::keep_interactions(complex, smallm, peptides, dist_cutoff);

        // Output phase
        lemon::print_residue_name_counts(std::cout, pdbid, complex, smallm);
    };

    current_path(p);
    lemon::run_archive(worker, vec, ncpu, chun);
}
