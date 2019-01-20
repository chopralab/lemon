from lemon import *

distance = 6.0

class MyWorkflow(Workflow):
    def worker(self, frame, pdbid):
        heme_names = ResidueNameSet()
        heme_names.append(ResidueName("HEM"))
        heme_names.append(ResidueName("HEA"))
        heme_names.append(ResidueName("HEB"))
        heme_names.append(ResidueName("HEC"))

        hemegs = select_specific_residues(frame, heme_names)
        smallm = select_small_molecules(frame, small_molecule_types, 10)

        # Pruning phase
        prune_identical_residues(frame, smallm)
        prune_cofactors(frame, smallm, common_cofactors)
        prune_cofactors(frame, smallm, common_fatty_acids)

        keep_interactions(frame, smallm, hemegs, distance)

        # Output phase
        return print_residue_name_counts(pdbid, frame, smallm)

    def finalize(self):
        pass
