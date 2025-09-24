#!/usr/bin/env python3
"""
Superposed Docking Analysis: WT vs Mutant Comparison

This script performs comparative analysis of molecular docking results from 
Wild-Type (WT) and Mutant protein systems. It aligns structures, identifies 
overlapping ligand poses, and generates visualization-ready PDB files.

Author: [Your Name]
Date: [Date]
Version: 1.0
"""

import mdtraj as md
import numpy as np
import glob
import os
from itertools import combinations


def load_structures(path_pattern, system_name="System"):
    """
    Load all PDB structures from specified path pattern.
    
    Args:
        path_pattern (str): Glob pattern for PDB files (e.g., "/path/*/*.pdb")
        system_name (str): Name for logging purposes
        
    Returns:
        list: List of (filename, trajectory) tuples
    """
    files = sorted(glob.glob(path_pattern))
    structures = []
    
    print(f"🔍 Loading {len(files)} structures from {system_name}...")
    
    for f in files:
        try:
            traj = md.load(f)
            structures.append((os.path.basename(f), traj))
        except Exception as e:
            print(f"⚠️ Error loading {f}: {e}")
            continue
    
    return structures


def get_alignment_indices(traj, selection="resname CU"):
    """
    Get atom indices for structural alignment.
    
    Args:
        traj (mdtraj.Trajectory): Input trajectory
        selection (str): MDTraj selection string for alignment atoms
        
    Returns:
        numpy.array: Array of atom indices
    """
    return traj.topology.select(selection)


def align_structures(structures, ref_traj, alignment_selection="resname CU"):
    """
    Align all structures to a reference structure.
    
    Args:
        structures (list): List of (name, trajectory) tuples
        ref_traj (mdtraj.Trajectory): Reference trajectory for alignment
        alignment_selection (str): Atom selection for alignment
        
    Returns:
        list: List of aligned (name, trajectory) tuples
    """
    ref_indices = get_alignment_indices(ref_traj, alignment_selection)
    
    if len(ref_indices) == 0:
        raise ValueError(f"❌ No atoms found for selection: {alignment_selection}")
    
    aligned_structures = []
    
    for name, traj in structures:
        indices = get_alignment_indices(traj, alignment_selection)
        
        if len(indices) != len(ref_indices):
            print(f"⚠️ {name}: Different number of alignment atoms "
                  f"({len(indices)} vs {len(ref_indices)}). Skipping.")
            continue
            
        try:
            aligned = traj.superpose(ref_traj, 
                                   atom_indices=indices, 
                                   ref_atom_indices=ref_indices)
            aligned_structures.append((name, aligned))
        except Exception as e:
            print(f"⚠️ Error aligning {name}: {e}")
            continue
    
    return aligned_structures


def extract_ligand_coords(traj, ligand_name="DMP"):
    """
    Extract ligand coordinates from trajectory.
    
    Args:
        traj (mdtraj.Trajectory): Input trajectory
        ligand_name (str): Residue name of the ligand
        
    Returns:
        numpy.array or None: Ligand coordinates in nm
    """
    ligand_indices = traj.topology.select(f"resname {ligand_name}")
    if len(ligand_indices) == 0:
        return None
    return traj.xyz[0, ligand_indices]


def calculate_ligand_rmsd(coords1, coords2):
    """
    Calculate RMSD between two sets of ligand coordinates.
    
    Args:
        coords1, coords2 (numpy.array): Coordinate arrays
        
    Returns:
        float: RMSD in Angstroms
    """
    if coords1 is None or coords2 is None:
        return float('inf')
    
    if len(coords1) != len(coords2):
        return float('inf')
    
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1))) * 10  # Convert nm to Å
    return rmsd


def write_combined_pdb(filename, protein_traj, wt_ligands, mut_ligands, ligand_name="DMP"):
    """
    Write combined PDB file with protein and all ligands.
    
    Args:
        filename (str): Output filename
        protein_traj (mdtraj.Trajectory): Protein structure
        wt_ligands (list): WT ligand structures
        mut_ligands (list): Mutant ligand structures
        ligand_name (str): Ligand residue name
    """
    with open(filename, 'w') as f:
        f.write("REMARK Combined protein structure with WT and Mutant ligands\n")
        f.write(f"REMARK WT ligands: {len(wt_ligands)} (chain W)\n")
        f.write(f"REMARK Mutant ligands: {len(mut_ligands)} (chain M)\n")
        
        atom_num = 1
        
        # Write protein atoms
        for i, atom in enumerate(protein_traj.topology.atoms):
            coord = protein_traj.xyz[0, i]
            chain_id = atom.residue.chain.chain_id or "A"
            element_symbol = (atom.element.symbol if atom.element is not None else " X").rjust(2)
            
            f.write(f"ATOM  {atom_num:5d} {atom.name:<4s} {atom.residue.name:>3s} {chain_id:1s}"
                   f"{atom.residue.resSeq:4d}    "
                   f"{coord[0]*10:8.3f}{coord[1]*10:8.3f}{coord[2]*10:8.3f}"
                   f"{1.00:6.2f}{20.00:6.2f}          {element_symbol}\n")
            atom_num += 1
        
        # Write WT ligands (chain W)
        for res_num, (name, traj) in enumerate(wt_ligands, 1):
            ligand_indices = traj.topology.select(f"resname {ligand_name}")
            if len(ligand_indices) > 0:
                for idx in ligand_indices:
                    atom = traj.topology.atom(idx)
                    coord = traj.xyz[0, idx]
                    element_symbol = (atom.element.symbol if atom.element is not None else " X").rjust(2)
                    
                    f.write(f"ATOM  {atom_num:5d} {atom.name:<4s} {ligand_name} W{res_num:4d}    "
                           f"{coord[0]*10:8.3f}{coord[1]*10:8.3f}{coord[2]*10:8.3f}"
                           f"{1.00:6.2f}{20.00:6.2f}          {element_symbol}\n")
                    atom_num += 1
                
                f.write(f"TER   {atom_num:5d}      {ligand_name} W{res_num:4d}\n")
                atom_num += 1
        
        # Write Mutant ligands (chain M)
        for res_num, (name, traj) in enumerate(mut_ligands, 1):
            ligand_indices = traj.topology.select(f"resname {ligand_name}")
            if len(ligand_indices) > 0:
                for idx in ligand_indices:
                    atom = traj.topology.atom(idx)
                    coord = traj.xyz[0, idx]
                    element_symbol = (atom.element.symbol if atom.element is not None else " X").rjust(2)
                    
                    f.write(f"ATOM  {atom_num:5d} {atom.name:<4s} {ligand_name} M{res_num:4d}    "
                           f"{coord[0]*10:8.3f}{coord[1]*10:8.3f}{coord[2]*10:8.3f}"
                           f"{1.00:6.2f}{20.00:6.2f}          {element_symbol}\n")
                    atom_num += 1
                
                f.write(f"TER   {atom_num:5d}      {ligand_name} M{res_num:4d}\n")
                atom_num += 1
        
        f.write("END\n")


def filter_overlapping_ligands(wt_aligned, mut_aligned, rmsd_threshold=3.0, ligand_name="DMP"):
    """
    Filter out overlapping ligands between WT and Mutant systems.
    
    Args:
        wt_aligned (list): Aligned WT structures
        mut_aligned (list): Aligned Mutant structures  
        rmsd_threshold (float): RMSD threshold in Angstroms
        ligand_name (str): Ligand residue name
        
    Returns:
        tuple: (wt_filtered, mut_filtered) - Non-overlapping structures
    """
    # Extract ligand coordinates
    wt_ligands = []
    mut_ligands = []
    
    for name, traj in wt_aligned:
        coords = extract_ligand_coords(traj, ligand_name)
        if coords is not None:
            wt_ligands.append((name, traj, coords))
    
    for name, traj in mut_aligned:
        coords = extract_ligand_coords(traj, ligand_name)
        if coords is not None:
            mut_ligands.append((name, traj, coords))
    
    print(f"🔍 Analyzing overlap between {len(wt_ligands)} WT and {len(mut_ligands)} Mutant ligands")
    
    # Find overlapping pairs
    overlapping_pairs = []
    overlapping_wt_indices = set()
    overlapping_mut_indices = set()
    
    for i, (wt_name, wt_traj, wt_coords) in enumerate(wt_ligands):
        for j, (mut_name, mut_traj, mut_coords) in enumerate(mut_ligands):
            rmsd = calculate_ligand_rmsd(wt_coords, mut_coords)
            
            if rmsd < rmsd_threshold:
                overlapping_pairs.append((i, j, rmsd, wt_name, mut_name))
                overlapping_wt_indices.add(i)
                overlapping_mut_indices.add(j)
                print(f"⚠️ Overlap detected: {wt_name} vs {mut_name} (RMSD: {rmsd:.2f} Å)")
    
    # Filter non-overlapping structures
    wt_filtered = [
        (wt_ligands[i][0], wt_ligands[i][1]) 
        for i in range(len(wt_ligands)) 
        if i not in overlapping_wt_indices
    ]
    
    mut_filtered = [
        (mut_ligands[j][0], mut_ligands[j][1]) 
        for j in range(len(mut_ligands)) 
        if j not in overlapping_mut_indices
    ]
    
    # Print summary
    print(f"✅ Filtering completed:")
    print(f"   - Original WT ligands: {len(wt_ligands)}")
    print(f"   - Non-overlapping WT ligands: {len(wt_filtered)}")
    print(f"   - WT ligands removed: {len(wt_ligands) - len(wt_filtered)}")
    print(f"   - Original Mutant ligands: {len(mut_ligands)}")
    print(f"   - Non-overlapping Mutant ligands: {len(mut_filtered)}")
    print(f"   - Mutant ligands removed: {len(mut_ligands) - len(mut_filtered)}")
    print(f"   - Total overlapping pairs: {len(overlapping_pairs)}")
    
    return wt_filtered, mut_filtered


def create_output_files(ref_traj, wt_aligned, mut_aligned, wt_filtered, mut_filtered,
                       output_prefix="docking_analysis", ligand_name="DMP"):
    """
    Create output PDB files.
    
    Args:
        ref_traj (mdtraj.Trajectory): Reference structure
        wt_aligned, mut_aligned (list): All aligned structures
        wt_filtered, mut_filtered (list): Filtered structures
        output_prefix (str): Prefix for output files
        ligand_name (str): Ligand residue name
    """
    # Extract protein structure
    protein_indices = ref_traj.topology.select("protein or resname CU")
    protein_traj = ref_traj.atom_slice(protein_indices)
    
    # Filter structures with ligands
    wt_with_ligands = [
        (name, traj) for name, traj in wt_aligned
        if len(traj.topology.select(f"resname {ligand_name}")) > 0
    ]
    
    mut_with_ligands = [
        (name, traj) for name, traj in mut_aligned
        if len(traj.topology.select(f"resname {ligand_name}")) > 0
    ]
    
    # Create output files
    all_ligands_file = f"{output_prefix}_all_ligands.pdb"
    filtered_file = f"{output_prefix}_filtered_ligands.pdb"
    
    # Complete dataset
    write_combined_pdb(all_ligands_file, protein_traj, wt_with_ligands, mut_with_ligands, ligand_name)
    print(f"✅ Complete dataset saved: {all_ligands_file}")
    print(f"   - WT ligands (chain W): {len(wt_with_ligands)}")
    print(f"   - Mutant ligands (chain M): {len(mut_with_ligands)}")
    
    # Filtered dataset
    write_combined_pdb(filtered_file, protein_traj, wt_filtered, mut_filtered, ligand_name)
    print(f"✅ Filtered dataset saved: {filtered_file}")
    print(f"   - Non-overlapping WT ligands: {len(wt_filtered)}")
    print(f"   - Non-overlapping Mutant ligands: {len(mut_filtered)}")


def main():
    """Main analysis workflow."""
    
    # Configuration
    config = {
        'wt_path_pattern': "/path/to/WT/*/*.pdb",  # Update this path
        'mut_path_pattern': "/path/to/MUTANT/*/*.pdb",  # Update this path
        'alignment_selection': "resname CU",
        'ligand_name': "DMP",
        'rmsd_threshold': 3.0,
        'output_prefix': "docking_analysis"
    }
    
    print("🚀 Starting Superposed Docking Analysis")
    print("=" * 50)
    
    # Step 1: Load structures
    print("📁 Step 1: Loading structures...")
    wt_structures = load_structures(config['wt_path_pattern'], "WT")
    mut_structures = load_structures(config['mut_path_pattern'], "Mutant")
    
    if not wt_structures or not mut_structures:
        raise ValueError("❌ No structures loaded. Check your path patterns.")
    
    print(f"✅ Loaded: {len(wt_structures)} WT | {len(mut_structures)} Mutant structures")
    
    # Step 2: Structural alignment
    print("\n🧬 Step 2: Performing structural alignment...")
    ref_name, ref_traj = wt_structures[0]
    print(f"🎯 Using '{ref_name}' as reference structure")
    
    wt_aligned = align_structures(wt_structures, ref_traj, config['alignment_selection'])
    mut_aligned = align_structures(mut_structures, ref_traj, config['alignment_selection'])
    
    print(f"✅ Alignment completed: {len(wt_aligned)} WT | {len(mut_aligned)} Mutant")
    
    # Step 3: Filter overlapping ligands
    print(f"\n🔍 Step 3: Filtering overlapping ligands (RMSD < {config['rmsd_threshold']} Å)...")
    wt_filtered, mut_filtered = filter_overlapping_ligands(
        wt_aligned, mut_aligned, 
        config['rmsd_threshold'], 
        config['ligand_name']
    )
    
    # Step 4: Generate output files
    print(f"\n📄 Step 4: Generating output files...")
    create_output_files(
        ref_traj, wt_aligned, mut_aligned, 
        wt_filtered, mut_filtered,
        config['output_prefix'], config['ligand_name']
    )
    
    print("\n🎉 Analysis completed successfully!")
    print("=" * 50)
    print("📋 Summary:")
    print(f"   - Total WT structures processed: {len(wt_aligned)}")
    print(f"   - Total Mutant structures processed: {len(mut_aligned)}")
    print(f"   - Non-overlapping WT ligands: {len(wt_filtered)}")
    print(f"   - Non-overlapping Mutant ligands: {len(mut_filtered)}")
    print(f"   - RMSD threshold used: {config['rmsd_threshold']} Å")


if __name__ == "__main__":
    # Example usage - update paths as needed
    try:
        main()
    except Exception as e:
        print(f"❌ Error: {e}")
        print("Please check your configuration and input files.")
