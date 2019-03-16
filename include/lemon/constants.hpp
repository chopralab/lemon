#ifndef LEMON_CONSTANTS_HPP
#define LEMON_CONSTANTS_HPP

#include <unordered_set>

#include "lemon/residue_name.hpp"

namespace lemon {

//! Common peptides in the PDB
//!
//! Names of the most common peptides found in the PDB. This list includes the
//! twenty natural amino acids as well as common post-translation modifications.
//! These were Identified using the 'peptide_residues' **Lemon** workflow.
//!
//! 'Non-natural' amino acids:
//! - *CSD*: 3-sulfinoalanine
//! - *PCA*: pyroglutamic acid
//! - *DLE*: D-leucine
//! - *KCX*: lysine nz-carboxylic acid
//! - *CAS*: S-(dimethylarsenic) cysteine
//! - *CSO*: S-hydroxycysteine
//! - *PTR*: O-phosphotyrosine
//! - *TPO*: phosphothreonine
//! - *SEP*: phosphoserine
//! - *CME*: S,S-(2-hydroxyethyl)thiocysteine
//! - *SAH*: S-adenosyl-l-homocysteine
//! - *MLY*: N-dimethyl-lysine
//! - *HYP*: 4-hydroxyproline
//! - *MSE*: selenomethionine
const ResidueNameSet common_peptides{
    {"CSD"}, {"PCA"}, {"DLE"}, {"KCX"}, {"CAS"}, {"CSO"}, {"PTR"},
    {"CME"}, {"SAH"}, {"TPO"}, {"SEP"}, {"MLY"}, {"HYP"}, {"MSE"},
    {"CYS"}, {"TRP"}, {"MET"}, {"HIS"}, {"TYR"}, {"GLN"}, {"PHE"},
    {"ASN"}, {"PRO"}, {"ARG"}, {"THR"}, {"ASP"}, {"ILE"}, {"LYS"},
    {"SER"}, {"GLU"}, {"VAL"}, {"GLY"}, {"ALA"}, {"LEU"},
};

//! Common co-factors in the PDB
//!
//! Many proteins require a cofactor to function properly. These cofators are
//! present with the corresponding protein. Common cofactors include: flavin
//! groups (*FAD*, *FMN*), nicotinamide groups (*NAD*, *NAP*), heme groups (
//! *HEM*, *HEA*, *HEB*, *HEC*), nucleotides (*ADP*, *ATP*, *GDP*, *GTP*),
//! citric acid (*CIT*, *FLC*), and S-adenosylmethionine (*SAM*).
const ResidueNameSet common_cofactors{
    {"FAD"}, {"FMN"}, {"NAD"}, {"NAP"}, {"CLA"}, {"HEM"}, {"HEA"}, {"HEB"},
    {"HEC"}, {"ADP"}, {"ATP"}, {"GDP"}, {"GTP"}, {"UNL"}, {"CIT"}, {"FLC"},
    {"BE7"}, {"MHA"}, {"DHD"}, {"B3P"}, {"BTB"}, {"NHE"}, {"GOL"}, {"DTP"},
    {"SAM"}, {"SIA"}, {"ICT"}, {"EPE"}, {"MES"},
};

//! Common fatty-acids in the PDB
//!
//! Several protein crystal structures contain fatty acids used to help induce
//! crystallization.  This list contains fatty acids and other 'linear' small
//! molecules.
const ResidueNameSet common_fatty_acids{
    {"PG6"}, {"PE7"}, {"PG5"}, {"PEU"}, {"PGE"}, {"PIG"}, {"PE8"}, {"PE4"},
    {"P33"}, {"C8E"}, {"OTE"}, {"XPE"}, {"N8E"}, {"DR6"}, {"PEG"}, {"2PE"},
    {"P6G"}, {"1PE"}, {"SPM"}, {"SPK"}, {"SPD"}, {"1PG"}, {"PG4"}, {"MYR"},
    {"OLA"}, {"OLB"}, {"OLC"}, {"PLM"}, {"PEE"}, {"LHG"}, {"MC3"}, {"PAM"},
};

const ResidueNameSet proline_res{
    {"PRO"},
    {"HYP"},
    {"PCA"},
};

//! Linkage types for small-molecules in the PDB
//!
//! This set contains the types of linkages used to initially identify a small
//! molecule in a crystal structure.*PEPTIDE-LIKE* is included here becasue
//! many drugs are classified with this linkage.
const std::unordered_set<std::string> small_molecule_types{
    {"NON-POLYMER"}, {"OTHER"}, {"PEPTIDE-LIKE"}};

//! PDB entries that should be skipped as they are time consuming to calculate.
//!
//! These PDB entries are of virus capsids and therefore are time consuming to
//! perform operations on.  The typical user may wish to skip these entries
//! as to reduce the total amount of time required for a **Lemon** workflow.
//! These entries are indentified as they take an order of magnitude longer to
//! perform an operation than other entries.
const std::unordered_set<std::string> large_entries{
    {"3J3Q"}, {"3J3Y"}, {"5Y6P"}};

} // namespace lemon

#endif
