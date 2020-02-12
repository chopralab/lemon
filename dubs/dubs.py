from candiy_lemon import lemon
import sys

# List of dictionaries to keep track of parts of the file

# Key: reference pdbid, Value: path to .mmtf file
pathDict = {}
# Key: reference pdbID, Value: chemical ID of ligand bound to reference protein for removal
referenceLigandDict = {}
# Key: reference pdbID, Value: list of sm or non-sm protein
referenceDict = {}
# Key: reference pdbID, Value: list of proteins to aling to reference (like in pinc)
alignProtDict = {}
# Key: align protein pdbID, Value: ligand chemical ID code for ligand removal 
alignProtLigandDict = {}
# Key: pdbID, Value: chemical id for SM ligand
pdbIDSMDict = {}
# Key: pdbID, Value: tuple(resCode, chainID, residue ID)
pdbIDNonSMDict = {}
# Key: pdbID, Value: chemical id for SM ligand
noAlignSMDict = {}
# Key: pdbID, Value: tuple(resCode, chainID, residue ID)
noAlignNonSMDict = {}

entries = lemon.Entries()

# Method for parsing a formated input file
def parse_input_file(fname):
    # Open file and initialize flags to 0
    f = open(fname,"r")
    curRefPdbID = ""
    flag = 0

    for line in f:
        # Check to see if the line contains any of the tags
        # Set appropriate flags if it does
        if line.startswith("@<reference>"):
            flag = 1
        elif line.startswith("@<align_prot>"):
            flag = 2
        elif line.startswith("@<align_sm_ligands>"):
            flag = 3
        elif line.startswith("@<align_non_sm_ligands>"):
            flag = 4
        elif line.startswith("@<no_align_sm_ligands>"):
            flag = 5
        elif line.startswith("@<no_align_non_sm_ligands>"):
            flag = 6
        elif line.startswith("@<end>"):
            flag = 0
        else:
            # If the line does not contain a flag
            # Add info to appropriate dictionary based of set flags
            if flag == 1:
                pdbID = line.split(" ")[0].strip().upper()
                path = line.split(" ")[1].strip()
                curRefPdbID = pdbID
                pathDict[pdbID] = path

                if len(line.split(" ")) == 3:
                    chemID = line.split(" ")[2].strip().upper()
                    referenceLigandDict[pdbID] = chemID 
            
            elif flag == 2:
                pdbID = line.split(" ")[0].strip()
                chemID = line.split(" ")[1].strip()

                if alignProtDict.get(curRefPdbID,0) == 0:
                    alignProtDict[curRefPdbID] = [pdbID]
                else:
                    alignProtDict[curRefPdbID].append(pdbID)

                if alignProtLigandDict.get(pdbID, 0) == 0:
                    alignProtLigandDict[pdbID] = [chemID]
                else:
                    alignProtLigandDict[pdbID].append(chemID)

                entries.add(pdbID)

            elif flag == 3:
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

            elif flag == 4:
                pdbID = line.split(" ")[0].strip().upper()
                residueCode = line.split(" ")[1].split("-")[0].strip().upper()
                chainID = line.split(" ")[1].split("-")[1].strip().upper()
                residueID = line.split(" ")[1].split("-")[2].strip().upper()

                if referenceDict.get(curRefPdbID,0) == 0:
                    referenceDict[curRefPdbID] = [pdbID]
                else:
                    referenceDict[curRefPdbID].append(pdbID)
                
                if pdbIDNonSMDict.get(pdbID,0) == 0:
                    pdbIDNonSMDict[pdbID] = [tuple([residueCode,chainID,residueID])]
                else:
                    pdbIDNonSMDict[pdbID].append(tuple([residueCode,chainID,residueID]))

                entries.add(pdbID)
            
            elif flag == 5:
                pdbID = line.split(" ")[0].strip().upper()
                chemID = line.split(" ")[1].strip().upper()

                if noAlignSMDict.get(pdbID,0) == 0:
                    noAlignSMDict[pdbID] = [chemID]
                else:
                    noAlignSMDict[pdbID].append(chemID)

                entries.add(pdbID)

            elif flag == 6:
                pdbID = line.split(" ")[0].strip().upper()
                residueCode = line.split(" ")[1].split("-")[0].strip().upper()
                chainID = line.split(" ")[1].split("-")[1].strip()
                residueID = line.split(" ")[1].split("-")[2].strip()

                if noAlignNonSMDict.get(pdbID,0) == 0:
                    noAlignNonSMDict[pdbID] = [tuple([residueCode,chainID,residueID])]
                else:
                    noAlignNonSMDict[pdbID].append(tuple([residueCode,chainID,residueID]))

                entries.add(pdbID)

# Get from the command line
# for testing we also can hard set the path
if len(sys.argv) > 4:
    input_file_path = sys.argv[1]
    hadoop_path = sys.argv[2]
    cores = int(sys.argv[3])
    outdir = sys.argv[4]
else:
    #TODO change this if needed for testing
    print("You must give the input file, path to RCSB hadoop files, number of cores, and working directory")
    sys.exit(1)

# Define Lemon workflow class
class MyWorkflow(lemon.Workflow):
    def __init__(self):
        lemon.Workflow.__init__(self)
        
        self.reference_structures = {}
        for key, value in pathDict.items():
            self.reference_structures[key] = lemon.open_file(value)

        self.noAlignSMDict = noAlignSMDict

        self.outdir = outdir

        self.write_all_proteins = False

    def worker(self, entry, pdbid):
        # Define and assign the reference pdbid
        refpdbid = ""
        # mode is 0 unassigned, 1 for alignment for protein, 2 for alignment for ligand
        mode = 0

        # Check for pdbID as a protein to be aligned (like in PINC)
        for key, value in alignProtDict.items():
            if pdbid in value:
                refpdbid = key
                mode = 1

        # Check for protein-ligand pair for alignment 
        for key, value in referenceDict.items():
            if pdbid in value:
                refpdbid = key
                mode = 2

        if mode == 0:
            for ligand_code in self.noAlignSMDict.get(pdbid, []):
                rns = lemon.ResidueNameSet()
                rns.append(lemon.ResidueName(ligand_code))
                ligand_ids = lemon.select_specific_residues(entry, rns)
                lemon.prune_identical_residues(entry, ligand_ids)

                for ligand_id in ligand_ids:
                    protein = lemon.Frame()
                    ligand = lemon.Frame()

                    lemon.separate_protein_and_ligand(entry, ligand_id, 25.0, protein, ligand)
                    lemon.write_file(protein, self.outdir + "/" + pdbid + "_" + ligand_code + ".pdb")
                    lemon.write_file(ligand, self.outdir + "/" + pdbid + "_" + ligand_code + ".sdf")

            return pdbid + " no alignment\n"

        elif mode == 1:
            # If we need to align to a protein (like in PINC)
            alignment = lemon.TMscore(entry, self.reference_structures[refpdbid])
            positions = entry.positions()
            lemon.align(positions, alignment.affine)

            lemon.write_file(entry, self.outdir + "/" + refpdbid + "_" + pdbid + ".pdb")

            return "Align Protein: " + refpdbid + "_" + pdbid + " to " + refpdbid + " with score of " + str(alignment.score) + "\n"

        elif mode == 2:
            # If we are doing ligand alignment, that can be done here
            # Get a list of the ligands associated with the protein we are trying to align
            SM_ligandList = pdbIDSMDict.get(pdbid, []) 
            Non_SM_ligandList = pdbIDNonSMDict.get(pdbid, [])

            alignment = lemon.TMscore(entry, self.reference_structures[refpdbid])
            positions = entry.positions()
            lemon.align(positions, alignment.affine)

            if len(SM_ligandList) > 0:

                for ligand_code in SM_ligandList:

                    rns = lemon.ResidueNameSet()
                    rns.append(lemon.ResidueName(ligand_code))

                    ligand_ids = lemon.select_specific_residues(entry, rns)
                    lemon.prune_identical_residues(entry, ligand_ids)

                    for ligand_id in ligand_ids:
                        protein = lemon.Frame()
                        ligand = lemon.Frame()

                        lemon.separate_protein_and_ligand(entry, ligand_id, 25.0, protein, ligand)
                        lemon.write_file(ligand, self.outdir + "/" + refpdbid + "_" + pdbid + "_" + ligand_code + ".sdf")

                        if self.write_all_proteins:
                            lemon.write_file(protein, self.outdir + "/" + refpdbid + "_" + pdbid + "_" + ligand_code + ".pdb")

            return "Align Protein: " + pdbid + " to " + refpdbid + " with score of " + str(alignment.score) + "\n"

    def finalize(self):
        pass

# Parse the input file
parse_input_file(input_file_path)
print(pathDict)
# Initilize the workflow
wf = MyWorkflow()

# TODO Get these from the command-line or ask the user
lemon.launch(wf, hadoop_path, cores, entries)
