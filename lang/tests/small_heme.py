from lemon import *

distance = 6.0

class MyWorkflow(Workflow):
    def worker(self, entry, pdbid):
        heme_names = ResidueNameSet()
        heme_names.append(ResidueName("HEM"))
        heme_names.append(ResidueName("HEA"))
        heme_names.append(ResidueName("HEB"))
        heme_names.append(ResidueName("HEC"))

        hemegs = select_specific_residues(entry, heme_names)
        smallm = select_small_molecules(entry, small_molecule_types, 10)

        # Pruning phase
        prune_identical_residues(entry, smallm)
        prune_cofactors(entry, smallm, common_cofactors)
        prune_cofactors(entry, smallm, common_fatty_acids)

        keep_interactions(entry, smallm, hemegs, distance)

        # Output phase
        return print_residue_name_counts(pdbid, entry, smallm)

    def finalize(self):
        pass
