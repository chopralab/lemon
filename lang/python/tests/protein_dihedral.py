from __future__ import print_function
import lemon

class MyWorkflow(lemon.Workflow):
    def __init__(self):
        import lemon
        lemon.Workflow.__init__(self)

        self.dihedral_dict = {}

    def worker(self, entry, pdbid):
        import lemon
        import math

        protein_only = lemon.Frame()
        peptides = lemon.ResidueIDs()

        if (lemon.select_specific_residues(entry, peptides,
                                           lemon.common_peptides) == 0):
            return ""

        lemon.separate_residues(entry, peptides, protein_only)
        dihedrals = protein_only.topology().dihedrals()

        for dihedral in dihedrals:
            dihedralnm = ""
            try:
                dihedralnm = lemon.protein_dihedral_name(protein_only, dihedral, lemon.proline_res)
            except lemon.GeometryError as error:
                return pdbid + ": " + 'error' + '\n'

            theta = protein_only.dihedral(dihedral[0], dihedral[1],
                                          dihedral[2], dihedral[3])

            dbin = int(math.floor(theta / 0.01))

            sbin = (dihedralnm, dbin)

            if sbin in self.dihedral_dict:
                self.dihedral_dict[sbin] = self.dihedral_dict[sbin] + 1
            else:
                self.dihedral_dict[sbin] = 1

        return ""

    def finalize(self):
        for sbin, count in self.dihedral_dict.items():
            print(sbin[0], '\t', sbin[1], '\t', count)

wf = MyWorkflow()

lemon.launch(wf, LEMON_HADOOP_DIR, LEMON_NUM_THREADS)
