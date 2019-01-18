from lemon import *

native = open_file("../../test/files/1AAQ.mmtf")

class MyWorkflow(Workflow):
    def worker(self, frame, pdbid):
        junk = PositionVec()

        tm = TMscore(frame, native, junk, False)

        return pdbid + "\t" + str(tm.score) + "\t" + str(tm.rmsd) + "\t" + str(tm.aligned) + "\n"
    def finalize(self):
        pass
