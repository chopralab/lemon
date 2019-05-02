#include <iostream>
#include <sstream>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"
#include "lemon/tmalign.hpp"

int main(int argc, char* argv[]) {
    lemon::Options o;
    auto distance = 2.4;
    o.add_option("--distance,-d", distance,
                 "Largest distance between a metal and a small molecule.");
    std::string reference("ref.pdb");
    o.add_option("--reference,-r", reference, "Reference DMFase");
    o.parse_command_line(argc, argv);

    auto ref = chemfiles::Trajectory(reference).read();
    auto ref_fe1 = lemon::select::specific_residues(ref, {"FE"});
    auto ref_res = lemon::select::specific_residues(ref, {"GLU", "TYR"});

    // Only the first iron from this file
    lemon::prune::keep_interactions(ref, ref_res, {*ref_fe1.begin()}, distance);

    std::vector<chemfiles::Vector3D> x;
    for (auto res : ref_res) {
        auto cur_res = ref.topology().residues()[res];
        auto loc = lemon::tmalign::find_element_by_name(ref, cur_res, "CA");
        x.push_back(ref.positions()[loc]);
    }

    std::vector<double> w(x.size(), 1.0);

    auto worker = [distance, &x, w](chemfiles::Frame entry,
                                    const std::string& pdbid) {

        // Selection phase
        auto metals = lemon::select::metal_ions(entry);
        //auto metals = lemon::select::specific_residues(entry, {"FE"});
        auto resids = lemon::select::specific_residues(entry,
            {"GLU", "ASP", "TYR"}
        );

        std::string result;

        if (metals.size() == 0) {
            return result;
        }

        for (auto metal : metals) {
            auto resids2 = resids;
            lemon::prune::keep_interactions(entry, resids2, {metal}, distance);
            std::vector<chemfiles::Vector3D> y;

            if (resids2.size() != 3) {
                continue;
            }

            result += entry.topology().residues()[metal].name() + " ";

            for (auto res : resids2) {
                auto cur_res = entry.topology().residues()[res];
                auto loc = lemon::tmalign::find_element_by_name(entry, cur_res, "CA");
                y.push_back(entry.positions()[loc]);
                result += cur_res.name() + std::to_string(*cur_res.id()) + ":" + cur_res.get("chainname")->as_string() + " ";
            }

            chemfiles::Matrix3D u = chemfiles::Matrix3D::unit();
            chemfiles::Vector3D t;
            int err;

            auto rms = lemon::tmalign::kabsch(w, x,
                {y[0], y[1], y[2]},
                x.size(), u, t, err
            );
            result += std::to_string(rms) + " ";

            rms = lemon::tmalign::kabsch(w, x,
                {y[0], y[2], y[1]},
                x.size(), u, t, err
            );
            result += std::to_string(rms) + " ";

            rms = lemon::tmalign::kabsch(w, x,
                {y[1], y[0], y[2]},
                x.size(), u, t, err
            );
            result += std::to_string(rms) + " ";

            rms = lemon::tmalign::kabsch(w, x,
                {y[1], y[2], y[0]},
                x.size(), u, t, err
            );
            result += std::to_string(rms) + " ";

            rms = lemon::tmalign::kabsch(w, x,
                {y[2], y[0], y[1]},
                x.size(), u, t, err
            );
            result += std::to_string(rms) + " ";

            rms = lemon::tmalign::kabsch(w, x,
                {y[2], y[1], y[0]},
                x.size(), u, t, err
            );
            result += std::to_string(rms) + " ";

            result += pdbid + "\n";
        }

        // Output phase
        return result;
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
