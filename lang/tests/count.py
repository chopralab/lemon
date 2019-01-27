from __future__ import print_function
from lemon import *

count = []

class MyWorkflow(Workflow):
    def worker(self, entry, pdbid):
        smallm = select_molecules(entry, small_molecule_types, 10)
        select_metal_ions(entry, smallm)
        count.append(smallm.size())
    def finalize(self):
        print(str(count))
