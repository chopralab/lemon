#ifndef SEPARATE_HPP
#define SEPARATE_HPP

#include "chemfiles/Frame.hpp"

namespace benchmarker {
    void separate_protein_and_ligand(
        const chemfiles::Frame& input,
        const chemfiles::Residue& ligand_residue,
        chemfiles::Frame& protein,
        chemfiles::Frame& ligand,
        double pocket_size = 0
    );
}

#endif
