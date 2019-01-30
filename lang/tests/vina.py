import lemon

class MyWorkflow(lemon.Workflow):
    def worker(self, entry, pdbid):
        import lemon

        # Selection phase
        smallm = lemon.select_small_molecules(entry, lemon.small_molecule_types, 10)
        if (len(smallm) == 0):
            return ""

        # Pruning phase
        lemon.prune_identical_residues(entry, smallm)
        lemon.prune_cofactors(entry, smallm, lemon.common_cofactors)
        lemon.prune_cofactors(entry, smallm, lemon.common_fatty_acids)

        # Output phase
        prot = lemon.ResidueIDs()
        residues = entry.topology().residues()
        for resid in range(0, len(residues)):
            prot.append(resid)
        
        result = ""
        for smallm_id in smallm:
            lig_copy = lemon.ResidueIDs()
            lig_copy.append(smallm_id)

            # Hack to remove self
            prot_copy = lemon.ResidueIDs(prot)
            lemon.remove_interactions(entry, prot_copy, lig_copy, 0.001)

            vscore = lemon.vina_score(entry, smallm_id, prot_copy, 8.0)

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
