from lemon import *

class MyWorkflow(Workflow):
    def worker(self, frame, pdbid):
        # Selection phase
        smallm = select_small_molecules(frame, small_molecule_types, 10)
        if (smallm.size() == 0):
            return ""

        # Pruning phase
        prune_identical_residues(frame, smallm)
        prune_cofactors(frame, smallm, common_cofactors)
        prune_cofactors(frame, smallm, common_fatty_acids)

        # Output phase
        prot = ResidueIDs()
        residues = frame.topology().residues()
        for resid in range(0, residues.size()):
            prot.append(resid)
        
        result = ""
        for smallm_id in smallm:
            lig_copy = ResidueIDs()
            lig_copy.append(smallm_id)

            # Hack to remove self
            prot_copy = ResidueIDs(prot)
            remove_interactions(frame, prot_copy, lig_copy, 0.001)

            vscore = vina_score(frame, smallm_id, prot_copy, 8.0)

            result += pdbid + "\t" 
            result += residues[smallm_id].name() + "\t" 
            result += str(vscore.g1) + "\t" 
            result += str(vscore.g2) + "\t"
            result += str(vscore.hydrogen) + "\t" 
            result += str(vscore.hydrophobic) + "\t"
            result += str(vscore.rep) + "\n"

        return result
    def finalize(self):
        pass
