from candiy_lemon import lemon
import sys, os
import numpy as np
import pandas as pd

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

# Maximum distance allowed between a residue and the ligand
maximumResidueDistance = 25.0

entries = set()

# Method for getting binding affinity information
# Input is list of tuple with each tuple -> (pdb_id,lig_code)
def get_bind_affinty(bind_tup_list):
    # We can get the csv file directly from a url using pandas
    # Use custom columns due to some of the columns not being needed
    binding_moad_url = "http://bindingmoad.org/files/csv/every_bind.csv"
    binding_moad_data = pd.read_csv(binding_moad_url,usecols=[0,2,3,4,5,7,8,9],header=None)

    # Set the columns manually because there is no header given
    binding_moad_data.columns = ["ec","pdb_id","lig_code","validity","affnty_type","affnty_val","affnty_units","SMILES"]

    # Create our dictionaries and varaibles here
    pdb_id_dict = {}
    lig_id_dict = {}
    cur_pdb_id = ""

    # Go through each of the row and add the data to the dictionaries
    for index,row in binding_moad_data.iterrows():
        if pd.isna(row["pdb_id"]) == False:
            cur_pdb_id = row["pdb_id"]
        elif (pd.isna(row["lig_code"]) == False) and (pd.isna(row["affnty_type"]) == False) and (pd.isna(row["affnty_val"]) == False) and (pd.isna(row["affnty_units"]) == False):
            simple_lig_code = row["lig_code"].split(":")[0]
            simple_lig_code = simple_lig_code.split(" ")

            for code in simple_lig_code:

                if pdb_id_dict.get(cur_pdb_id,0) == 0:
                    pdb_id_dict[cur_pdb_id] = [code]
                else:
                    if simple_lig_code not in pdb_id_dict[cur_pdb_id]:
                        pdb_id_dict[cur_pdb_id].append(code)

                if lig_id_dict.get((cur_pdb_id,code),0) == 0:
                    lig_id_dict[(cur_pdb_id,code)] = [[row["lig_code"],row["affnty_type"],row["affnty_val"],row["affnty_units"]]]
                else:
                    lig_id_dict[(cur_pdb_id,code)].append([row["lig_code"],row["affnty_type"],row["affnty_val"],row["affnty_units"]])

    #To access the ligand we need to first check to see if the ligand exists in pdb_id_dict
    # If this is one of the ligands availibe, then we access the tuple in the lig_id_dict with the binding affinity
    ret_dict = {}
    for pair in bind_tup_list:
        pdb_id = pair[0]
        lig_code = pair[1]

        if lig_id_dict.get((pdb_id,lig_code),0) == 0:
            ret_dict[(pdb_id,lig_code)] = "No Binding Information Found!"
        else:
            ret_dict[(pdb_id,lig_code)] = lig_id_dict[(pdb_id,lig_code)]

    # Return the final dict, key -> (pdb_id,lig_code) and value is information on that ligand interaction
    return ret_dict

# Method for parsing a formated input file
def parse_input_file(fname):
    # Open file and initialize flags to 0
    f = open(fname,"r")
    curRefPdbID = ""
    flag = 0

    global maximumResidueDistance
    global referenceLigandDict
    global referenceDict
    global alignProtDict

    for line in f:
        # Check to see if the line contains any of the tags
        # Set appropriate flags if it does
        if line.startswith("#"):
            continue
        elif line.startswith("@<reference>"):
            flag = 1
            continue
        elif line.startswith("@<align_prot>"):
            flag = 2
            continue
        elif line.startswith("@<align_sm_ligands>"):
            flag = 3
            continue
        elif line.startswith("@<align_non_sm_ligands>"):
            flag = 4
            continue
        elif line.startswith("@<no_align_sm_ligands>"):
            flag = 5
            continue
        elif line.startswith("@<no_align_non_sm_ligands>"):
            flag = 6
            continue
        elif line.startswith("@<maximum_residue_distance>"):
            flag = 7
            continue
        elif line.startswith("@<end>"):
            flag = 0
            continue
        elif line.startswith("@"):
            print("Invalid tag (@<>) detected, please check all tags are in the proper format")
            sys.exit()

        # If the line does not contain a flag
        # Add info to appropriate dictionary based of set flags
        if flag == 1:
            # Error Check
            if len(line.split(" ")) < 2:
                print("Invalid input under @<reference> tag, please check to ensure proper input")
                sys.exit(1)

            pdbID = line.split(" ")[0].strip().upper()
            path = line.split(" ")[1].strip()
            curRefPdbID = pdbID
            pathDict[pdbID] = path

            if len(line.split(" ")) >= 3:
                chemID = line.split(" ")[2].strip().upper()
                referenceLigandDict[pdbID] = chemID 
            
        elif flag == 2:
            if len(line.split(" ")) < 2:
                print("Invalid input under @<align_prot> tag, please check to ensure proper input")
                sys.exit(1)

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
            if len(line.split(" ")) < 2:
                print("Invalid input under @<align_sm_ligands> tag, please check to ensure proper input")
                sys.exit(1)

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
            if len(line.split(" ")) < 3:
                print("Invalid input under @<align_non_sm_ligands> tag, please check to ensure proper input")
                sys.exit(1)

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
            if len(line.split(" ")) < 2:
                print("Invalid input under @<no_align_sm_ligands> tag, please check to ensure proper input")
                sys.exit(1)

            pdbID = line.split(" ")[0].strip().upper()
            chemID = line.split(" ")[1].strip().upper()

            if noAlignSMDict.get(pdbID,0) == 0:
                noAlignSMDict[pdbID] = [chemID]
            else:
                noAlignSMDict[pdbID].append(chemID)

            entries.add(pdbID)

        elif flag == 6:
            if len(line.split(" ")) < 3:
                print("Invalid input under @<no_align_non_sm_ligands> tag, please check to ensure proper input")
                sys.exit(1)

            pdbID = line.split(" ")[0].strip().upper()
            residueCode = line.split(" ")[1].split("-")[0].strip().upper()
            chainID = line.split(" ")[1].split("-")[1].strip()
            residueID = line.split(" ")[1].split("-")[2].strip()

            if noAlignNonSMDict.get(pdbID,0) == 0:
                noAlignNonSMDict[pdbID] = [tuple([residueCode,chainID,residueID])]
            else:
                noAlignNonSMDict[pdbID].append(tuple([residueCode,chainID,residueID]))

            entries.add(pdbID)
        
        elif flag == 7:
            maximumResidueDistance = float(line.strip())

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

        self.maximumResidueDistance = maximumResidueDistance

    def extact_ligand_write(self, entry, pdbid, ligand_code):
        ligand_name = ligand_code
        ligand_chain = ""
        ligand_altloc = "A"

        if ligand_code.find(":") != -1:
            ligand_split = ligand_code.split(":")
            ligand_chain = ligand_split[0]
            ligand_name = ligand_split[1]

        if ligand_name.find("_") != -1:
            ligand_split = ligand_name.split("_")
            ligand_name = ligand_split[0]
            ligand_altloc = ligand_split[1]

        ligand_ids = ""

        if ligand_name.startswith("+"):
            ids_to_use = set()
            for my_id in ligand_name.split("+"):
                if my_id.isdigit():
                    ids_to_use.add(int(my_id))
 
            ligand_ids = lemon.select_residue_ids(entry, ids_to_use)

            ligand_ids = lemon.prune_identical_residues(entry, ligand_ids)
            if ligand_chain:
                ligand_ids = lemon.has_property(entry, ligand_ids, "chainname", lemon.Property(ligand_chain))
                    
            protein = lemon.Frame()
            ligand = lemon.Frame()

            if not os.path.isdir(self.outdir + "/" + pdbid):
                os.mkdir(self.outdir + "/" + pdbid)

            lemon.separate_protein_and_ligands(entry, ligand_ids, self.maximumResidueDistance, protein, ligand, ligand_altloc)
            out_base = self.outdir + "/" + pdbid + "/" + pdbid + "_" + ligand_name

            lemon.write_file(protein, out_base + ".pdb")
            lemon.write_file(ligand, out_base + ".sdf")

            return

        rns = set()
        rns.add(lemon.ResidueName(ligand_name))
        ligand_ids = lemon.select_specific_residues(entry, rns)

        if ligand_chain:
            ligand_ids = lemon.has_property(entry, ligand_ids, "chainname", lemon.Property(ligand_chain))

        ligand_ids = lemon.prune_identical_residues(entry, ligand_ids)

        for ligand_id in ligand_ids:
            protein = lemon.Frame()
            ligand = lemon.Frame()

            if not os.path.isdir(self.outdir + "/" + pdbid):
                os.mkdir(self.outdir + "/" + pdbid)

            lemon.separate_protein_and_ligand(entry, ligand_id, self.maximumResidueDistance, protein, ligand, ligand_altloc)
            out_base = self.outdir + "/" + pdbid + "/" + pdbid + "_" + ligand_name

            lemon.write_file(protein, out_base + ".pdb")
            lemon.write_file(ligand, out_base + ".sdf")

            if len(ligand_ids) > 1:
                print("WARNING:It appears that {} in {}".format(ligand_code, pdbid), file=sys.stderr)
                print("        is defined multiple times in the same biological assembly,", file=sys.stderr)
                print("        only keeping the first time it appears", file=sys.stderr)
                break

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
                self.extact_ligand_write(entry, pdbid, ligand_code)
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
                    self.extact_ligand_write(entry, pdbid, ligand_code)
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
