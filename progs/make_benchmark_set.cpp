#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <chemfiles.hpp>

#include "lemon/entries.hpp"
#include "lemon/count.hpp"
#include "lemon/prune.hpp"
#include "lemon/run.hpp"
#include "lemon/select.hpp"
#include "lemon/separate.hpp"

using namespace boost::filesystem;

int main(int argc, char* argv[]) {
    path entries(argc > 1 ? argv[1] : "entries.idx");
    path p(argc > 2 ? argv[2] : ".");
    path o(argc > 3 ? absolute(argv[3]) : absolute("."));
    double dist_cutoff = argc > 4 ? std::strtod(argv[4], nullptr) : 10.0;
    size_t ncpu = argc > 5 ? std::strtoul(argv[5], nullptr, 0) : 1;
    size_t chun = argc > 6 ? std::strtoul(argv[6], nullptr, 0) : 1;

    if (!is_regular_file(entries)) {
        std::cerr << "You must supply a valid entries file" << std::endl;
        return 1;
    }

    if (!is_directory(p)) {
        std::cerr << "You must supply a valid directory" << std::endl;
        return 2;
    }

    lemon::PDBIDVec vec;
    std::unordered_map<std::string, lemon::ResidueNameSet> rnms;
    std::ifstream is(entries.string());
    lemon::read_entry_file(is, vec, rnms);

    auto worker = [dist_cutoff, &rnms, &o](const chemfiles::Frame& complex,
                                           const std::string& pdbid) {

        // Selection phase
        auto smallm = lemon::select_specific_residues(complex, rnms[pdbid]);

        // Pruning phase
        lemon::remove_identical_residues(complex, smallm);

        // Output phase
        for (auto resid : smallm) {
            chemfiles::Frame prot;
            chemfiles::Frame lig;
            lemon::separate_protein_and_ligand(complex, resid, prot, lig,
                                               dist_cutoff);

            path protfile = o;
            protfile /= pdbid + "_" + lig.get("name")->as_string() + ".pdb";
            path ligfile  = o;
            ligfile /= pdbid + "_" + lig.get("name")->as_string() + ".sdf";

            chemfiles::Trajectory prot_traj(protfile.string(),'w');
            chemfiles::Trajectory lig_traj(ligfile.string(),'w');
            prot_traj.write(prot);
            lig_traj.write(lig);
            prot_traj.close();
            lig_traj.close();
        }
    };

    current_path(p);
    lemon::call_multithreaded(worker, vec, ncpu, chun);
}
