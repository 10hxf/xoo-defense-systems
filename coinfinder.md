# Coinfinder analysis of co-occurrence and avoidance relationships

working path: /root/phagework/hxf-206xoo/coinfinder/

## Step 1: **Pan-genome analysis using Roary and phylogenetic tree construction using FastTree**

```
roary -e --mafft -p 24 gff/*.gff -f 206xooroary
```

Construct a phylogenetic tree using FastTree based on the Roary-generated `core_gene_alignment.aln` file

```
FastTree –nt –gtr core_gene_alignment.aln > 206xoo.newick
```

## Step 2: Prepare the defense system presence/absence list

The first column contains the defense system name, and the second column contains the strain ID, which corresponds to the filename of the GFF input file used when running Roary

Run the `extract_defense_systems.py` script to generate the `defense_systems_presence.tsv` file from the PADLOC results

```
#!/usr/bin/env python3
"""
Script to extract defense systems from PADLOC results and create a 
tab-delimited file for Coinfinder analysis.

This script processes PADLOC CSV output files for multiple strains and 
creates a two-column TSV file with defense system names and strain IDs.
"""

import os
import pandas as pd
import glob
from pathlib import Path

def extract_defense_systems_from_csv(csv_file_path):
    """
    Extract unique defense systems from a PADLOC CSV file.
    
    Args:
        csv_file_path (str): Path to the PADLOC CSV file
        
    Returns:
        set: Set of unique defense system names
    """
    try:
        # Read the CSV file
        df = pd.read_csv(csv_file_path)
        
        # Check if the required columns exist
        if 'system' not in df.columns:
            print(f"Warning: 'system' column not found in {csv_file_path}")
            return set()
        
        # Extract unique defense systems
        # Group by system.number to identify complete systems
        if 'system.number' in df.columns:
            # Get unique systems by grouping by system.number
            unique_systems = df.groupby('system.number')['system'].first().unique()
        else:
            # Fallback: just get unique system names
            unique_systems = df['system'].unique()
        
        # Remove any NaN values and convert to set
        defense_systems = set(system for system in unique_systems if pd.notna(system))
        
        return defense_systems
        
    except Exception as e:
        print(f"Error processing {csv_file_path}: {str(e)}")
        return set()

def process_all_strains(padloc_dir, output_file):
    """
    Process all strain directories and create the defense systems presence file.
    
    Args:
        padloc_dir (str): Path to the directory containing all strain subdirectories
        output_file (str): Path to the output TSV file
    """
    
    # List to store all defense system entries
    defense_entries = []
    
    # Get all strain directories
    strain_dirs = [d for d in os.listdir(padloc_dir) 
                   if os.path.isdir(os.path.join(padloc_dir, d))]
    
    print(f"Found {len(strain_dirs)} strain directories")
    
    processed_strains = 0
    total_systems = 0
    
    # Process each strain directory
    for strain_id in sorted(strain_dirs):
        strain_path = os.path.join(padloc_dir, strain_id)
        
        # Find CSV files in the strain directory
        csv_files = glob.glob(os.path.join(strain_path, "*.csv"))
        
        if not csv_files:
            print(f"Warning: No CSV files found in {strain_path}")
            continue
        
        # Process all CSV files for this strain (usually just one)
        strain_defense_systems = set()
        
        for csv_file in csv_files:
            defense_systems = extract_defense_systems_from_csv(csv_file)
            strain_defense_systems.update(defense_systems)
        
        # Add entries for this strain
        for defense_system in strain_defense_systems:
            defense_entries.append([defense_system, strain_id])
        
        if strain_defense_systems:
            processed_strains += 1
            total_systems += len(strain_defense_systems)
            print(f"Processed {strain_id}: {len(strain_defense_systems)} defense systems")
        else:
            print(f"No defense systems found for {strain_id}")
    
    # Create DataFrame and save to file
    if defense_entries:
        df_output = pd.DataFrame(defense_entries, columns=['Defense_System', 'Strain_ID'])
        
        # Sort by defense system name, then by strain ID
        df_output = df_output.sort_values(['Defense_System', 'Strain_ID'])
        
        # Save as tab-delimited file without header (as required by Coinfinder)
        df_output.to_csv(output_file, sep='\t', header=False, index=False)
        
        print(f"\nSummary:")
        print(f"- Processed {processed_strains} strains")
        print(f"- Total defense system instances: {len(defense_entries)}")
        print(f"- Unique defense systems: {len(df_output['Defense_System'].unique())}")
        print(f"- Output saved to: {output_file}")
        
        # Display some statistics
        print(f"\nDefense system frequency:")
        system_counts = df_output['Defense_System'].value_counts()
        print(system_counts.head(10))
        
    else:
        print("No defense systems found in any strain!")

def main():
    """Main function to run the defense systems extraction."""
    
    # Define paths
    padloc_directory = "/root/phagework/hxf-206xoo/coinfinder/padloc/"
    output_file = "defense_systems_presence.tsv"
    
    # Check if the input directory exists
    if not os.path.exists(padloc_directory):
        print(f"Error: Directory {padloc_directory} does not exist!")
        return
    
    print("Starting defense systems extraction...")
    print(f"Input directory: {padloc_directory}")
    print(f"Output file: {output_file}")
    print("-" * 60)
    
    # Process all strains
    process_all_strains(padloc_directory, output_file)
    
    print("-" * 60)
    print("Defense systems extraction completed!")

if __name__ == "__main__":
    main()
```

## Step 3: Co-occurrence and mutual exclusivity analysis

```
coinfinder -i defense_systems_presence.tsv -p 206xoo.newick -o defense_cooccurrence --associate -m -x 8
```

