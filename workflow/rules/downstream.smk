# rules/downstream.smk

# Conditional flags from config
run_mlst      = config["pipeline"].get("run_mlst", False)
run_amr       = config["pipeline"].get("run_amr", False)
run_virulence = config["pipeline"].get("run_virulence", False)
run_plasmid   = config["pipeline"].get("run_plasmid", False)

def get_contigs_for_downstream(wildcards):
    """Get contigs file for downstream analysis - either assembly output or outgroup FASTA"""
    outgroup_name = config["annotation"].get("outgroup_name", "outgroup")
    if wildcards.sample == outgroup_name:
        return config["annotation"]["outgroup_fasta"]
    else:
        return f"results/assembly/{wildcards.sample}/contigs.fasta"

# MLST Typing
rule mlst_typing:
    conda: 
        "../envs/downstream.yaml"
    threads: 2
    input:
        contigs = get_contigs_for_downstream
    output:
        "results/downstream/mlst/{sample}.txt"
    params:
        scheme = config["mlst"]["scheme"]
    log:
        "logs/mlst/{sample}.log"
    run:
        if not run_mlst:
            shell("echo 'Skipping MLST for {wildcards.sample}'")
        else:
            shell("mlst --scheme {params.scheme} {input.contigs} > {output} 2> {log}")

# AMR Screen
rule amr_screen:
    conda: 
        "../envs/downstream.yaml"
    threads: 4
    input:
        contigs = get_contigs_for_downstream
    output:
        "results/downstream/amr/{sample}.tsv"
    params:
        extra = config["amr"]["extra"]
    log:
        "logs/amr/{sample}.log"
    run:
        if not run_amr:
            shell("echo 'Skipping AMR for {wildcards.sample}'")
        else:
            shell("abricate --db card {params.extra} {input.contigs} > {output} 2> {log}")

# Virulence Screen
rule virulence_screen:
    conda: 
        "../envs/downstream.yaml"
    threads: 4
    input:
        contigs = get_contigs_for_downstream
    output:
        "results/downstream/virulence/{sample}.tsv"
    params:
        extra = config["virulence"]["extra"]
    log:
        "logs/virulence/{sample}.log"
    run:
        if not run_virulence:
            shell("echo 'Skipping virulence for {wildcards.sample}'")
        else:
            shell("abricate --db vfdb {params.extra} {input.contigs} > {output} 2> {log}")

# PLSDB Sketch
rule mash_sketch_plasmid_db:
    conda: 
        "../envs/downstream.yaml"
    input:
        fasta = config["plasmid_db"]["fasta"]
    output:
        "databases/plsdb/plasmids.msh"
    threads: 1
    shell:
        "mkdir -p databases/plsdb && mash sketch -o databases/plsdb/plasmids {input.fasta}"

# PLSDB BLAST DB
rule makeblastdb_plasmid:
    conda: 
        "../envs/downstream.yaml"
    input:
        fasta = config["plasmid_db"]["fasta"]
    output:
        nsq= "databases/plsdb/plasmids.nsq",
        nin= "databases/plsdb/plasmids.nin",
        nhr= "databases/plsdb/plasmids.nhr"
    shell:
        "mkdir -p databases/plsdb && makeblastdb -in {input.fasta} -dbtype nucl -out databases/plsdb/plasmids"

# Plasmid Mash Screen
rule mash_screen_plasmids:
    conda: 
        "../envs/downstream.yaml"
    threads: 4
    input:
        contigs = get_contigs_for_downstream,
        db = "databases/plsdb/plasmids.msh"
    output:
        "results/downstream/plasmid/{sample}/mash_hits.txt"
    params:
        thresh = config["plasmid_db"]["mash_threshold"]
    log:
        "logs/plasmid/{sample}.mash.log"
    run:
        if not run_plasmid:
            shell("touch {output}")
        else:
            shell("mash screen -p {threads} -i {params.thresh} {input.db} {input.contigs} > {output} 2> {log}")

# Plasmid BLASTn (mash-guided or direct)
rule plasmid_blastn:
    conda: 
        "../envs/downstream.yaml"
    threads: 4
    input:
        contigs    = get_contigs_for_downstream,
        hits       = "results/downstream/plasmid/{sample}/mash_hits.txt",
        nsq        = "databases/plsdb/plasmids.nsq"
    output:
        "results/downstream/plasmid/{sample}/blast_plasmid_hits.tsv"
    params:
        blast_db = "databases/plsdb/plasmids",
        evalue = config["plasmid_db"]["blast_evalue"],
        pid    = config["plasmid_db"]["blast_pid"]
    log:
        "logs/plasmid/{sample}.blast.log"
    run:
        if not run_plasmid:
            shell("touch {output}")
        else:
            # Check if mash found any hits - use size-based approach for all samples
            shell("""
                echo "Using size-based BLAST approach for plasmid detection" >> {log}
                
                # Strategy 1: Look for contigs in typical plasmid size range (1kb-500kb)
                python3 -c "
import sys
with open('{input.contigs}', 'r') as f:
    content = f.read()
    sequences = content.split('>')[1:]  # Skip empty first element
    candidates = []
    for seq in sequences:
        lines = seq.strip().split('\\n')
        header = lines[0]
        sequence = ''.join(lines[1:])
        length = len(sequence)
        if 1000 <= length <= 500000:  # Plasmid size range
            candidates.append(header.split()[0])  # Take first part of header
    
    with open('results/downstream/plasmid/{wildcards.sample}/candidates.txt', 'w') as out:
        for candidate in candidates[:15]:  # Limit to 15 candidates
            out.write(candidate + '\\n')
"
                
                # If no plasmid-sized contigs found, try smaller range
                if [ ! -s results/downstream/plasmid/{wildcards.sample}/candidates.txt ]; then
                    echo "No contigs in 1kb-500kb range, trying 500bp-50kb range" >> {log}
                    python3 -c "
import sys
with open('{input.contigs}', 'r') as f:
    content = f.read()
    sequences = content.split('>')[1:]
    candidates = []
    for seq in sequences:
        lines = seq.strip().split('\\n')
        header = lines[0]
        sequence = ''.join(lines[1:])
        length = len(sequence)
        if 500 <= length <= 50000:
            candidates.append(header.split()[0])
    
    with open('results/downstream/plasmid/{wildcards.sample}/candidates.txt', 'w') as out:
        for candidate in candidates[:20]:
            out.write(candidate + '\\n')
"
                fi
                
                # If still no candidates, try any contigs
                if [ ! -s results/downstream/plasmid/{wildcards.sample}/candidates.txt ]; then
                    echo "No suitable sized contigs, trying first 10 contigs" >> {log}
                    grep "^>" {input.contigs} | sed 's/^>//' | cut -d' ' -f1 | head -10 > results/downstream/plasmid/{wildcards.sample}/candidates.txt
                fi
                
                # Proceed with BLAST if we have candidates
                if [ -s results/downstream/plasmid/{wildcards.sample}/candidates.txt ]; then
                    echo "Found $(wc -l < results/downstream/plasmid/{wildcards.sample}/candidates.txt) candidates for BLAST" >> {log}
                    seqtk subseq {input.contigs} results/downstream/plasmid/{wildcards.sample}/candidates.txt > results/downstream/plasmid/{wildcards.sample}/query_seqs.fasta
                    
                    if [ -s results/downstream/plasmid/{wildcards.sample}/query_seqs.fasta ]; then
                        blastn -query results/downstream/plasmid/{wildcards.sample}/query_seqs.fasta \
                            -db {params.blast_db} -evalue {params.evalue} -perc_identity {params.pid} \
                            -outfmt '6 qseqid sseqid pident length evalue bitscore' -num_threads {threads} \
                            -max_target_seqs 2 > results/downstream/plasmid/{wildcards.sample}/blast_raw.tsv 2>> {log}
                        
                        # Filter for minimum alignment length (500bp) and sort by bitscore
                        awk '$4 >= 500 {{print}}' results/downstream/plasmid/{wildcards.sample}/blast_raw.tsv | \
                        sort -k6,6nr > {output}
                        
                        echo "BLAST completed, found $(wc -l < {output}) high-quality hits (after filtering)" >> {log}
                    else
                        echo "No sequences extracted for BLAST" >> {log}
                        touch {output}
                    fi
                else
                    echo "No candidates found for BLAST" >> {log}
                    touch {output}
                fi
            """)

# Aggregation Targets Helper
def define_agg_targets():
    targets = []
    if run_mlst:
        targets += expand("results/downstream/mlst/{sample}.txt", sample=DS_SAMPLES)
    if run_amr:
        targets += expand("results/downstream/amr/{sample}.tsv", sample=DS_SAMPLES)
    if run_virulence:
        targets += expand("results/downstream/virulence/{sample}.tsv", sample=DS_SAMPLES)
    if run_plasmid:
        targets += expand("results/downstream/plasmid/{sample}/blast_plasmid_hits.tsv", sample=DS_SAMPLES)
    return targets

# Downstream Aggregation Rule
rule aggregate_downstream:
    conda: 
        "../envs/downstream.yaml"
    input:
        entries = define_agg_targets()
    output:
        "results/downstream/summary.xlsx"
    script:
        "../scripts/aggregate_downstream.py"