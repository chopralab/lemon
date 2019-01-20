from __future__ import print_function
from lemon import *

rnc = ResidueNameCount()

class MyWorkflow(Workflow):
    def worker(self, frame, pdbid):
        count_residues(frame, rnc)
        return ""
    def finalize(self):
        for rn in rnc:
            print(str(rn.first) + '\t' + str(rn.second))
