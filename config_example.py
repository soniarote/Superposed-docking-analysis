"""
Configuration file for Superposed Docking Analysis

Copy this file and modify the paths and parameters according to your system.
"""

# Input paths - Update these to match your directory structure
WT_PATH_PATTERN = "/path/to/WT/*/*.pdb"
MUTANT_PATH_PATTERN = "/path/to/MUTANT/*/*.pdb"

# Analysis parameters
ALIGNMENT_SELECTION = "resname CU"  # Atoms used for structural alignment
LIGAND_NAME = "DMP"                 # Ligand residue name
RMSD_THRESHOLD = 3.0               # RMSD threshold in Angstroms for overlap detection

# Output settings  
OUTPUT_PREFIX = "docking_analysis"  # Prefix for output files

# Example configurations for different systems:

# Configuration 1: Metalloenzymes with copper binding sites
METALLOENZYME_CONFIG = {
    'wt_path_pattern': "/data/metalloenzyme/WT/docking_results/*/*.pdb",
    'mut_path_pattern': "/data/metalloenzyme/MUTANT/docking_results/*/*.pdb", 
    'alignment_selection': "resname CU",
    'ligand_name': "DMP",
    'rmsd_threshold': 3.0,
    'output_prefix': "metalloenzyme_analysis"
}

# Configuration 2: Standard protein-ligand systems
STANDARD_CONFIG = {
    'wt_path_pattern': "/data/protein/WT/poses/*/*.pdb",
    'mut_path_pattern': "/data/protein/MUT/poses/*/*.pdb",
    'alignment_selection': "protein and backbone",
    'ligand_name': "LIG", 
    'rmsd_threshold': 2.5,
    'output_prefix': "protein_ligand_analysis"
}

# Configuration 3: Enzyme active site focused
ACTIVE_SITE_CONFIG = {
    'wt_path_pattern': "/data/enzyme/WT/binding/*/*.pdb",
    'mut_path_pattern': "/data/enzyme/MUT/binding/*/*.pdb", 
    'alignment_selection': "resid 100 to 150 and name CA",  # Active site residues
    'ligand_name': "SUB",
    'rmsd_threshold': 2.0,
    'output_prefix': "active_site_analysis"
}
