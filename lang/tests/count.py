from __future__ import print_function
import lemon

class MyWorkflow(lemon.Workflow):
    def __init__(self):
        import lemon
        lemon.Workflow.__init__(self)
        self.count = []
    def worker(self, entry, pdbid):
        import lemon
        smallm = lemon.select_small_molecules(entry, lemon.small_molecule_types, 10)
        self.count.append(smallm.__len__())
        return ""
    def finalize(self):
        print(str(self.count))
