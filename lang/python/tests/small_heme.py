import lemon

class MyWorkflow(lemon.Workflow):
    def worker(self, entry, pdbid):
        import lemon
        heme_names = set()
        heme_names.add(lemon.ResidueName("HEM"))
        heme_names.add(lemon.ResidueName("HEA"))
        heme_names.add(lemon.ResidueName("HEB"))
        heme_names.add(lemon.ResidueName("HEC"))

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

wf = MyWorkflow()

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)
