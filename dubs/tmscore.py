from candiy_lemon import lemon

# List of dictionaries to keep track of parts of the file
pathDict = {}
referenceDict = {}
pdbIDSMDict = {}
pdbIDNonSMDict = {}

entries_to_use = lemon.Entries()

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

                entries.add(pdbID)

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

                entries.add(pdbID)

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
        # TODO load the reference files and use them in the worker thread
    def worker(self, entry, pdbid):
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

        tm = lemon.TMscore(entry, self.native, junk, False)
        print(pdbid + "\t" + str(tm.score) + "\t" + str(tm.rmsd) + "\t" + str(tm.aligned) + "\n")
        return pdbid + "\t" + str(tm.score) + "\t" + str(tm.rmsd) + "\t" + str(tm.aligned) + "\n"
    def finalize(self):
        pass


wf = MyWorkflow()

# TODO Get these from the command-line or ask the user
lemon.launch(wf, "../../full", 8)
