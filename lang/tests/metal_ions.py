from lemon import *

distance = 6.0

class MyWorkflow(Workflow):
    def worker(self, frame, pdbid):
        # Selection phase
        metals = select_metal_ions(frame)

        # Output phase
        return print_residue_name_counts(pdbid, frame, metals)
    def finalize(self):
        pass
