from lemon import *

class MyWorkflow(Workflow):
    def worker(self, entry, pdbid):
        import lemon
        heme_names = lemon.ResidueNameSet()
        heme_names.append(lemon.ResidueName("HEM"))
        heme_names.append(lemon.ResidueName("HEA"))
        heme_names.append(lemon.ResidueName("HEB"))
        heme_names.append(lemon.ResidueName("HEC"))

        hemegs = lemon.select_specific_residues(entry, heme_names)
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
