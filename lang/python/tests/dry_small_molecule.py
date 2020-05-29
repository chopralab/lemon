import lemon

class MyWorkflow(lemon.Workflow):
    def worker(self, entry, pdbid):
        import lemon
        wat_name = set()
        wat_name.add(lemon.ResidueName("HOH"))

        hemegs = lemon.select_specific_residues(entry, wat_name)
        smallm = lemon.select_small_molecules(entry, lemon.small_molecule_types, 10)

        # Pruning phase
        smallm = lemon.prune_identical_residues(entry, smallm)
        smallm = lemon.prune_cofactors(entry, smallm, lemon.common_cofactors)
        smallm = lemon.prune_cofactors(entry, smallm, lemon.common_fatty_acids)

        smallm = lemon.keep_interactions(entry, smallm, hemegs, 6.0)

        # Output phase
        return pdbid + lemon.count_print_residue_names(entry, smallm) + '\n'

    def finalize(self):
        pass

wf = MyWorkflow()

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)
