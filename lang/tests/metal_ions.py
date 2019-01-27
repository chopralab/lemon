from lemon import *

distance = 6.0

class MyWorkflow(Workflow):
    def worker(self, entry, pdbid):
        # Selection phase
        metals = select_metal_ions(entry)

        # Output phase
        return print_residue_name_counts(pdbid, entry, metals)
    def finalize(self):
        pass
