#!/usr/bin/env python3
"""
Prepare Phylogeny Input Script
Converts Roary core genome alignment to formats suitable for phylogenetic analysis
"""

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os

def convert_fasta_to_phylip(fasta_file, phylip_file):
    """Convert FASTA alignment to PHYLIP format for IQ-TREE"""
    try:
        # Read the alignment
        alignment = AlignIO.read(fasta_file, "fasta")
        
        # Write in PHYLIP format
        AlignIO.write(alignment, phylip_file, "phylip-relaxed")
        
        return True
    except Exception as e:
        print(f"Error converting to PHYLIP: {e}")
        return False

def analyze_alignment(fasta_file, info_file):
    """Analyze alignment and write statistics"""
    try:
        alignment = AlignIO.read(fasta_file, "fasta")
        
        n_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        # Calculate some basic statistics
        stats = {
            'sequences': n_sequences,
            'length': alignment_length,
            'total_sites': n_sequences * alignment_length
        }
        
        # Count variable sites (simplified analysis)
        variable_sites = 0
        for i in range(alignment_length):
            column = alignment[:, i]
            unique_chars = set(column.upper())
            # Remove gaps and ambiguous characters for counting
            unique_chars = {c for c in unique_chars if c in 'ATGC'}
            if len(unique_chars) > 1:
                variable_sites += 1
        
        stats['variable_sites'] = variable_sites
        stats['conserved_sites'] = alignment_length - variable_sites
        
        # Write information file
        with open(info_file, 'w') as f:
            f.write("Core Genome Alignment Statistics\n")
            f.write("=" * 40 + "\n")
            f.write(f"Number of sequences: {stats['sequences']}\n")
            f.write(f"Alignment length: {stats['length']:,} bp\n")
            f.write(f"Variable sites: {stats['variable_sites']:,}\n")
            f.write(f"Conserved sites: {stats['conserved_sites']:,}\n")
            f.write(f"Proportion variable: {stats['variable_sites']/stats['length']:.4f}\n")
            f.write("\nAlignment ready for phylogenetic analysis with IQ-TREE\n")
            f.write("Recommended IQ-TREE command:\n")
            f.write(f"iqtree -s {os.path.basename(info_file.replace('_info.txt', '.phy'))} -m MFP -bb 1000 -nt AUTO\n")
        
        return stats
        
    except Exception as e:
        print(f"Error analyzing alignment: {e}")
        return None

def main():
    # Get input and output files from Snakemake
    core_alignment = snakemake.input.core_alignment
    phylip_output = snakemake.output.phylip_alignment
    info_output = snakemake.output.alignment_info
    
    print(f"Converting core genome alignment to PHYLIP format...")
    print(f"Input: {core_alignment}")
    print(f"Output: {phylip_output}")
    
    # Convert to PHYLIP format
    success = convert_fasta_to_phylip(core_alignment, phylip_output)
    
    if success:
        print("Successfully converted to PHYLIP format")
        
        # Analyze alignment
        print("Analyzing alignment statistics...")
        stats = analyze_alignment(core_alignment, info_output)
        
        if stats:
            print(f"Alignment contains {stats['sequences']} sequences")
            print(f"Alignment length: {stats['length']:,} bp")
            print(f"Variable sites: {stats['variable_sites']:,} ({stats['variable_sites']/stats['length']:.2%})")
            print("Alignment is ready for phylogenetic analysis!")
        else:
            print("Warning: Could not analyze alignment statistics")
    else:
        print("Error: Failed to convert alignment")
        exit(1)

if __name__ == "__main__":
    main()
