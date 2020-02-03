import lemon

distance = 6.0

class MyWorkflow(lemon.Workflow):
    def worker(self, entry, pdbid):
        import lemon
        # Selection phase
        metals = lemon.select_metal_ions(entry)

        # Output phase
        return pdbid + lemon.count_print_residue_names(entry, metals) + '\n'
    def finalize(self):
        pass

wf = MyWorkflow()

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)
