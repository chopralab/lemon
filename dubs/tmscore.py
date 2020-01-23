from candiy_lemon import lemon

# Define constants used later in Lemon test script
LEMON_HADOOP_DIR = "../../full"
LEMON_NUM_THREADS = 8

# List of dictonaries to keep track of parts of the file
pathDict = {}
referenceDict = {}
pdbIDSMDict = {}
pdbIDNonSMDict = {}

# Method for parsing a formated input file
def parse_input_file(fname):
    f = open(fname,"r")
    curRefPdbID = ""
    refFlag = 0
    SMLigFlag = 0
    nonSMligFlag = 0
    for line in f:
        if line.startswith("@<reference>"):
            refFlag = 1
        elif line.startswith("@<align_sm_ligands>"):
            refFlag = 0
            SMLigFlag = 1
        elif line.startswith("@<align_non_sm_ligands>"):
            SMLigFlag = 0
            nonSMligFlag = 1
        elif line.startswith("@<end>"):
            nonSMligFlag = 0
        else:
            if refFlag == 1:
                pdbID = line.split(" ")[0].strip()
                path = line.split(" ")[1].strip()
                curRefPdbID = pdbID
                pathDict[pdbID] = path
            elif SMLigFlag == 1:
                pdbID = line.split(" ")[0].strip()
                chemID = line.split(" ")[1].strip()

                if referenceDict.get(curRefPdbID,0) == 0:
                    referenceDict[curRefPdbID] = [pdbID]
                else:
                    referenceDict[curRefPdbID].append(pdbID)
                
                if pdbIDSMDict.get(pdbID,0) == 0:
                    pdbIDSMDict[pdbID] = [chemID]
                else:
                    pdbIDSMDict[pdbID].append(chemID)

            elif nonSMligFlag == 1:
                pdbID = line.split(" ")[0].strip()
                chainID = line.split(" ")[1].split("-")[0].strip()
                residNum = line.split(" ")[1].split("-")[1].strip()

                if referenceDict.get(curRefPdbID,0) == 0:
                    referenceDict[curRefPdbID] = [pdbID]
                else:
                    referenceDict[curRefPdbID].append(pdbID)
                
                if pdbIDNonSMDict.get(pdbID,0) == 0:
                    pdbIDNonSMDict[pdbID] = [tuple([chainID,residNum])]
                else:
                    pdbIDNonSMDict[pdbID].append(tuple([chainID,residNum]))

"""
fname = "pinc_input.txt"
parse_input_file(fname)

print("Path Dict")
print(pathDict)
print("Reference Dict")
print(referenceDict)
print("SM Dict")
print(pdbIDSMDict)
print("Non SM Dict")
print(pdbIDNonSMDict)
"""

# Define Lemon workflow class
class MyWorkflow(lemon.Workflow):
    def __init__(self):
        lemon.Workflow.__init__(self)
        #self.native = lemon.open_file("../../../test/files/1AAQ.mmtf")
    def worker(self, entry, pdbid):
        # If the pdbid is not used, then we just skip over it
        if pdbid not in pdbIDSMDict.keys() and pdbid not in pdbIDNonSMDict.keys():
            return ""
        
        junk = lemon.PositionVec()

        # Define and assign the reference pdbid
        refpdbid = ""
        for key, value in referenceDict.items():
            if pdbid in value:
                refpdbid = key

        # Get the path to the reference file and set native to it
        refPath = pathDict[refpdbid]
        self.native = lemon.open_file(refPath)

        # Get a list of the ligands associated with the protein we are trying to align
        SM_ligandList = pdbIDSMDict[pdbid] 
        Non_SM_ligandList = pdbIDNonSMDict[pdbid]

        tm = le'mon.TMscore(entry, self.native, junk, False)
        print(pdbid + "\t" + str(tm.score) + "\t" + str(tm.rmsd) + "\t" + str(tm.aligned) + "\n")
        return pdbid + "\t" + str(tm.score) + "\t" + str(tm.rmsd) + "\t" + str(tm.aligned) + "\n"
    def finalize(self):
        pass


wf = MyWorkflow()

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)