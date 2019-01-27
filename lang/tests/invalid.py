from lemon import *

class MyWorkflow(Workflow):
    def worker(self, entry, pdbid, junk, junk2):
        smallm = select_molecules(entry, small_molecule_types, 10)
        select_metal_ions(entry, smallm)
        grr.append(smallm.size())
