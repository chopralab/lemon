from lemon import *

distance = 6.0

class MyWorkflow(Workflow):
    def worker(self, entry, pdbid):
        smallm = select_small_molecules(entry, small_molecule_types, 10)

        # Pruning phase
        prune_identical_residues(entry, smallm)
        prune_cofactors(entry, smallm, common_cofactors)
        prune_cofactors(entry, smallm, common_fatty_acids)

        # Output phase
        return print_residue_name_counts(pdbid, entry, smallm)
