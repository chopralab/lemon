#include <map>

#include <algorithm>
#include <iostream>
#include <sstream>

#include "parser.hpp"

#include "lemon/launch.hpp"
#include "lemon/lemon.hpp"
#include "lemon/tmalign.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o;
    std::string outdir = ".";
    auto distance = 6.0;
    auto dump_input = false;
    auto all_proteins = false;

    o.add_option("--distance,-d", distance,
                 "Largest distance between protein and a small molecule.");
    o.add_option("--outdir,-o", outdir, "output directory");
    o.add_flag("--all_proteins", all_proteins, "Include all proteins, even for "
        "complexes tagged as small molecule only");
    o.add_flag("--dump_input", dump_input, "Dump the input file (for debugging)");
    o.parse_command_line(argc, argv);

    outdir += "/";

    std::ifstream is;
    is.open(o.entries(), std::ios::in);

    DUBSParser parser;
    try {
        parser.parse(is);
    } catch (const std::exception& e) {
        std::cerr << "Could not parse " << o.entries() << ": " << e.what() << std::endl;
        return 1;
    }

    if (dump_input) {
        std::cout << parser.dump() << std::endl;
        return 0;
    }

    auto worker = [distance, all_proteins, &outdir, &parser]
        (chemfiles::Frame entry, const std::string& PDBid) {

        auto& rns = parser.ligands(PDBid);

        // Selection phase
        auto ligand_ids = lemon::select::specific_residues(entry, rns);

        if (parser.tag_type(PDBid) == DUBSParser::TAG_TYPE::REFERENCE) { // reference file
            auto protfile = outdir + PDBid + "_REF.pdb";

            chemfiles::Trajectory reference(protfile, 'w');
            reference.write(entry);

            return "Added reference file: " + PDBid + "\n";
        }

        // Pruning phase
        lemon::prune::identical_residues(entry, ligand_ids);

        auto reference = parser.reference(PDBid);
        auto result_str = PDBid;
        auto outdir_local = outdir;

        if (!reference.empty()) {
            auto& reference_struct = parser.reference_structure(reference);

            auto alignment = lemon::tmalign::TMscore(entry, reference_struct);

            auto pos = entry.positions();
            lemon::align(pos, alignment.affine);

            result_str += " aligned to " + reference + " with score of " +
                std::to_string(alignment.score) +
                "(" + std::to_string(alignment.aligned) + ")";

            auto name = parser.name(reference);
            if (!name.empty()) {
                outdir_local += name + "/" + name + "_";
            }
        }

        // Output phase
        for (auto resid : ligand_ids) {
            chemfiles::Frame prot;
            chemfiles::Frame lig;
            lemon::separate::protein_and_ligand(entry, resid, distance, prot,
                                                lig);

            auto lig_name = lig.get<chemfiles::Property::STRING>("name").value_or("UNK");

            result_str += " and ligand " + lig_name;

            auto lig_file = outdir_local + PDBid + "_" + lig_name + ".sdf";
            chemfiles::Trajectory lig_trj(lig_file, 'w');
            lig_trj.write(lig);

            if (!all_proteins && !reference.empty()) {
                continue;
            }

            auto prot_file = outdir_local + PDBid + "_" + lig_name + ".pdb";
            chemfiles::Trajectory prot_trj(prot_file, 'w');
            prot_trj.write(prot);
        }

        return result_str + "\n";
    };

    parser.make_directories(outdir);

    auto collector = lemon::print_combine(std::cout);

    try {
        lemon::run_parallel(worker, o.work_dir(), collector, o.ncpu(), parser.entries());
    } catch (std::runtime_error& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    return 0;
}
