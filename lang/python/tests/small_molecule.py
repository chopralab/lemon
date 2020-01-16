import lemon

class MyWorkflow(lemon.Workflow):
    def worker(self, entry, pdbid):
        import lemon
        smallm = lemon.select_small_molecules(entry, lemon.small_molecule_types, 10)

        # Pruning phase
        lemon.prune_identical_residues(entry, smallm)
        lemon.prune_cofactors(entry, smallm, lemon.common_cofactors)
        lemon.prune_cofactors(entry, smallm, lemon.common_fatty_acids)

        # Output phase
        return pdbid + lemon.count_print_residue_names(entry, smallm)

wf = MyWorkflow()

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)
