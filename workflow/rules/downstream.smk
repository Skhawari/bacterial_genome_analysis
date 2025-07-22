# rules/downstream.smk

# Conditional flags from config
run_mlst      = config["pipeline"].get("run_mlst", False)
run_amr       = config["pipeline"].get("run_amr", False)
run_virulence = config["pipeline"].get("run_virulence", False)
run_plasmid   = config["pipeline"].get("run_plasmid", False)

# MLST Typing
rule mlst_typing:
    conda: 
        "../envs/downstream.yaml"
    threads: 2
    input:
        contigs = "results/assembly/{sample}/contigs.fasta"
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
        contigs = "results/assembly/{sample}/contigs.fasta"
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
        contigs = "results/assembly/{sample}/contigs.fasta"
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
rule plasmid_mash_screen:
    conda: 
        "../envs/downstream.yaml"
    threads: 4
    input:
        contigs = "results/assembly/{sample}/contigs.fasta",
        db = "databases/plsdb/plasmids.msh"
    output:
        "results/downstream/plasmid/{sample}/mash_hits.txt"
    params:
        thresh = config["plasmid_db"]["mash_threshold"]
    log:
        "logs/plasmid/{sample}.mash.log"
    run:
        if not run_plasmid:
            shell("echo 'Skipping plasmid mash for {wildcards.sample}'")
        else:
            shell("mash screen -p {threads} -i {params.thresh} {input.db} {input.contigs} > {output} 2> {log}")

# Plasmid BLASTn
rule plasmid_blastn:
    conda: 
        "../envs/downstream.yaml"
    threads: 4
    input:
        contigs    = "results/assembly/{sample}/contigs.fasta",
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
            shell("echo 'Skipping plasmid blast for {wildcards.sample}'")
        else:
            shell("""
                awk '$1!~/^#/ {print $2}' {input.hits} | sort -u > downstream/plasmid/{wildcards.sample}/candidates.txt && \
                blastn -query <(seqtk subseq {input.contigs} downstream/plasmid/{wildcards.sample}/candidates.txt) \
                    -db {params.blast_db} -evalue {params.evalue} -perc_identity {params.pid} \
                    -outfmt '6 qseqid sseqid pident length evalue bitscore' -num_threads {threads} \
                    > {output} 2> {log}
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
        targets += expand("resutls/downstream/plasmid/{sample}/blast_plasmid_hits.tsv", sample=DS_SAMPLES)
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
        "scripts/aggregate_downstream.py"