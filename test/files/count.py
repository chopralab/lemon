from __future__ import print_function
from candiy_lemon import lemon

class MyWorkflow(lemon.Workflow):
    def __init__(self):
        lemon.Workflow.__init__(self)
        self.count = []
    def worker(self, entry, pdbid):
        smallm = lemon.select_small_molecules(entry, lemon.small_molecule_types, 10)
        self.count.append(smallm.__len__())
        return ""
    def finalize(self):
        print(str(self.count))

a = MyWorkflow()

lemon.launch(a, "rcsb_hadoop", 2)



