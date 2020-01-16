from __future__ import print_function
import lemon

class MyWorkflow(lemon.Workflow):
    def __init__(self):
        import lemon
        lemon.Workflow.__init__(self)
        self.rnc = lemon.ResidueNameCount()
    def worker(self, entry, pdbid):
        import lemon
        lemon.count_residues(entry, self.rnc)
        return ""
    def finalize(self):
        for rn in self.rnc:
            print(str(rn) + '\t' + str(self.rnc[rn]))

wf = MyWorkflow()

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)
