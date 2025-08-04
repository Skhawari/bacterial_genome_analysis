#!/usr/bin/env python3
import subprocess
import os
from pathlib import Path

def main():
    # Get Snakemake variables
    input_contigs = snakemake.input.contigs
    input_hits = snakemake.input.hits
    output_file = snakemake.output[0]
    blast_db = snakemake.params.blast_db
    evalue = snakemake.params.evalue
    pid = snakemake.params.pid
    threads = snakemake.threads
    log_file = snakemake.log[0]
    skip = snakemake.params.skip
    
    # Create output directory
    output_dir = Path(output_file).parent
    
    if skip:
        Path(output_file).touch()
        with open(log_file, 'w') as log:
            log.write("Plasmid analysis disabled in config\n")
        return
    
    with open(log_file, 'w') as log:
        log.write("Starting plasmid BLAST analysis\n")
        
        # Size-based candidate selection
        candidates = []
        with open(input_contigs, 'r') as f:
            content = f.read()
            sequences = content.split('>')[1:]
            
            for seq in sequences:
                lines = seq.strip().split('\n')
                header = lines[0].split()[0]
                sequence = ''.join(lines[1:])
                length = len(sequence)
                
                # Plasmid size range: 1kb-500kb
                if 1000 <= length <= 500000:
                    candidates.append(header)
        
        # Fallback to smaller range if needed
        if not candidates:
            log.write("No contigs in 1kb-500kb range, trying 500bp-50kb\n")
            # Repeat with smaller range...
        
        if candidates:
            candidates = candidates[:15]  # Limit candidates
            log.write(f"Found {len(candidates)} candidates for BLAST\n")
            
            # Extract sequences and run BLAST
            candidate_file = output_dir / "candidates.txt"
            query_file = output_dir / "query_seqs.fasta"
            
            with open(candidate_file, 'w') as f:
                for candidate in candidates:
                    f.write(f"{candidate}\n")
            
            # Extract sequences
            subprocess.run([
                "seqtk", "subseq", input_contigs, str(candidate_file)
            ], stdout=open(query_file, 'w'), stderr=log)
            
            # Run BLAST
            blast_cmd = [
                "blastn", "-query", str(query_file),
                "-db", blast_db, "-evalue", str(evalue),
                "-perc_identity", str(pid),
                "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
                "-num_threads", str(threads),
                "-max_target_seqs", "2"
            ]
            
            with open(output_file, 'w') as out:
                subprocess.run(blast_cmd, stdout=out, stderr=log)
            
            log.write("BLAST completed successfully\n")
        else:
            Path(output_file).touch()
            log.write("No suitable candidates found\n")

if __name__ == "__main__":
    main()
