from lemon import *

grr = []

class MyWorkflow(Workflow):
    def worker(self, frame):
        smallm = select_molecules(frame, small_molecule_types, 10)
        select_metal_ions(frame, smallm)
        grr.append(smallm.size())
    def finalize(self):
        print str(grr)
