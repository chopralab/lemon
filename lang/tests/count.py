from lemon import *

count = []

class MyWorkflow(Workflow):
    def worker(self, frame, pdbid):
        smallm = select_molecules(frame, small_molecule_types, 10)
        select_metal_ions(frame, smallm)
        count.append(smallm.size())
    def finalize(self):
        print str(count)
