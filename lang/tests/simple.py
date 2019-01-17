from lemon import *
class MyWorkflow(Workflow):
    def worker(self, frame):
        return get(frame.topology().residue(1), "chainname").get().as_string() + '\n'
    def finalize(self):
        pass
