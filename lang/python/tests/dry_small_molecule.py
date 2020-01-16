from lemon import *

class MyWorkflow(Workflow):
    def worker(self, entry, pdbid):
        import lemon
        wat_name = lemon.ResidueNameSet()
        wat_name.append(lemon.ResidueName("HOH"))

        hemegs = lemon.select_specific_residues(entry, wat_name)
        smallm = lemon.select_small_molecules(entry, lemon.small_molecule_types, 10)

        # Pruning phase
        lemon.prune_identical_residues(entry, smallm)
        lemon.prune_cofactors(entry, smallm, lemon.common_cofactors)
        lemon.prune_cofactors(entry, smallm, lemon.common_fatty_acids)

        lemon.keep_interactions(entry, smallm, hemegs, 6.0)

        # Output phase
        return pdbid + lemon.count_print_residue_names(entry, smallm) + '\n'

    def finalize(self):
        pass
