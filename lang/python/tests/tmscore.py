import lemon

class MyWorkflow(lemon.Workflow):
    def __init__(self):
        import lemon
        lemon.Workflow.__init__(self)
        self.native = lemon.open_file("../../../test/files/1AAQ.mmtf")
    def worker(self, entry, pdbid):
        import lemon

        tm = lemon.TMscore(entry, self.native)

        return pdbid + "\t" + str(tm.score) + "\t" + str(tm.aligned) + "\n"
    def finalize(self):
        pass

wf = MyWorkflow()

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)
