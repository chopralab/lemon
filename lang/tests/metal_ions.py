from lemon import *

distance = 6.0

class MyWorkflow(Workflow):
    def worker(self, entry, pdbid):
        import lemon
        # Selection phase
        metals = lemon.select_metal_ions(entry)

        # Output phase
        return pdbid + lemon.count_print_residue_names(entry, metals) + '\n'
    def finalize(self):
        pass
