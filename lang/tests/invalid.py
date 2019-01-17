from lemon import *

class MyWorkflow(Workflow):
    def worker(self, frame, pdbid, junk, junk2):
        smallm = select_molecules(frame, small_molecule_types, 10)
        select_metal_ions(frame, smallm)
        grr.append(smallm.size())
