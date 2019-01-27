from lemon import *
class MyWorkflow(Workflow):
    def worker(self, entry, pdbid):
        return entry.topology().residue(1).get("chainname").get().as_string() + '\n'
    def finalize(self):
        pass
