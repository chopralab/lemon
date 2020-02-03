from __future__ import print_function
import lemon
class MyWorkflow(lemon.Workflow):
    def worker(self, entry, pdbid):
        return '\t' + pdbid + '\n'
    def finalize(self):
        pass

wf = MyWorkflow()

print("Default arguments: ")

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)

# Re-run the work flow with a limited set of entries

limited = lemon.Entries()
limited.add("1DZF")
limited.add("1DZI")

print ("With selected Entries: " + str(limited))

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS, limited)
