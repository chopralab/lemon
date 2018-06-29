#include "benchmarker/parse.hpp"

namespace benchmarker {

void retreive_residue_counts(const chemfiles::Frame& file, ResnCount& resn_count) {
    auto& residues = file.topology().residues();

    for (auto& residue : residues) {
        auto resn_count_iterator = resn_count.find(residue.name());

        if (resn_count_iterator == resn_count.end()) {
            resn_count[residue.name()] = 1;
            continue;
        }

        ++resn_count_iterator->second;
    }
}

}
