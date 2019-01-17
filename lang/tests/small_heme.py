from lemon import *

distance = 6.0

class MyWorkflow(Workflow):
    def worker(self, frame, pdbid):
        smallm = select_small_molecules(frame, small_molecule_types, 10)

        # Pruning phase
        prune_identical_residues(frame, smallm)
        prune_cofactors(frame, smallm, common_cofactors)
        prune_cofactors(frame, smallm, common_fatty_acids)

        # Output phase
        return print_residue_name_counts(pdbid, frame, smallm)
