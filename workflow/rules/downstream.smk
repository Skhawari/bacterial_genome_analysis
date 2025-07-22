# rules/downstream.smk

run_mlst      = config["pipeline"].get("run_mlst", False)
run_amr       = config["pipeline"].get("run_amr", False)
run_virulence = config["pipeline"].get("run_virulence", False)
run_plasmid   = config["pipeline"].get("run_plasmid", False)

# MLST
rule mlst_typing:
    conda: "envs/downstream.yaml"
    threads: 2
    input:
        contigs = "results/assembly/{sample}/contigs.fasta"
    output:
        txt = "results/downstream/mlst/{sample}.txt"
    params:
        scheme = config["mlst"]["scheme"]
    log:
        "logs/mlst/{sample}.log"
    shell:
        """
        {{"echo 'Skipping MLST for {wildcards.sample}'" if not run_mlst else
        "mlst --scheme {params.scheme} {input.contigs} > {output.txt} 2> {log}"}}
        """

# AMR
rule amr_screen:
    conda: "envs/downstream.yaml"
    threads: 4
    input:
        contigs = "results/assembly/{sample}/contigs.fasta"
    output:
        tsv = "results/downstream/amr/{sample}.tsv"
    params:
        extra = config["amr"]["extra"]
    log:
        "logs/amr/{sample}.log"
    shell:
        """
        {{"echo 'Skipping AMR for {wildcards.sample}'" if not run_amr else
        "abricate --db card {params.extra} {input.contigs} > {output.tsv} 2> {log}"}}
        """

# Virulence
rule virulence_screen:
    conda: "envs/downstream.yaml"
    threads: 4
    input:
        contigs = "results/assembly/{sample}/contigs.fasta"
    output:
        tsv = "results/downstream/virulence/{sample}.tsv"
    params:
        extra = config["virulence"]["extra"]
    log:
        "logs/virulence/{sample}.log"
    shell:
        """
        {{"echo 'Skipping virulence for {wildcards.sample}'" if not run_virulence else
        "abricate --db vfdb {params.extra} {input.contigs} > {output.tsv} 2> {log}"}}
        """

# Mash Sketch
rule mash_sketch_plasmid_db:
    conda: "envs/downstream.yaml"
    input:
        fasta = config["plasmid_db"]["fasta"]
    output:
        msh = "databases/plsdb/plasmids.msh"
    threads: 1
    shell:
        "mash sketch -o databases/plsdb/plasmids {input.fasta}"

# BLAST DB
rule makeblastdb_plasmid:
    conda: "envs/downstream.yaml"
    input:
        fasta = config["plasmid_db"]["fasta"]
    output:
        touch("databases/plsdb/plasmids.nsq")
    shell:
        "makeblastdb -in {input.fasta} -dbtype nucl -out databases/plsdb/plasmids"

# Plasmid Mash
rule plasmid_mash_screen:
    conda: "envs/downstream.yaml"
    threads: 4
    input:
        contigs = "results/assembly/{sample}/contigs.fasta",
        db = "databases/plsdb/plasmids.msh"
    output:
        hits = "results/downstream/plasmid/{sample}/mash_hits.txt"
    params:
        thresh = config["plasmid_db"]["mash_threshold"]
    log:
        "logs/plasmid/{sample}.mash.log"
    shell:
        """
        {{"echo 'Skipping plasmid mash for {wildcards.sample}'" if not run_plasmid else
        "mash screen -p {threads} -i {params.thresh} {input.db} {input.contigs} > {output.hits} 2> {log}"}}
        """

# Plasmid BLASTn
rule plasmid_blastn:
    conda: "envs/downstream.yaml"
    threads: 4
    input:
        contigs = "results/assembly/{sample}/contigs.fasta",
        hits = "results/downstream/plasmid/{sample}/mash_hits.txt",
        blast_db = "databases/plsdb/plasmids"
    output:
        blast_result = "results/downstream/plasmid/{sample}/blast_plasmid_hits.tsv"
    params:
        evalue = config["plasmid_db"]["blast_evalue"],
        pid = config["plasmid_db"]["blast_pid"]
    log:
        "logs/plasmid/{sample}.blast.log"
    shell:
        """
        {{"echo 'Skipping plasmid blast for {wildcards.sample}'" if not run_plasmid else
        "awk '$1!~/^#/ {print $2}' {input.hits} | sort -u > results/downstream/plasmid/{wildcards.sample}/candidates.txt && "
        "seqtk subseq {input.contigs} results/downstream/plasmid/{wildcards.sample}/candidates.txt > results/downstream/plasmid/{wildcards.sample}/candidates.fa && "
        "blastn -query results/downstream/plasmid/{wildcards.sample}/candidates.fa "
        "-db {input.blast_db} -evalue {params.evalue} -perc_identity {params.pid} "
        "-outfmt '6 qseqid sseqid pident length evalue bitscore' -num_threads {threads} "
        "> {output.blast_result} 2> {log}"}}
        """

rule aggregate_downstream:
    conda: "../envs/downstream.yaml"
    input:
        entries = define_agg_targets()
    output:
        "results/downstream/summary.xlsx"
    script:
        "scripts/aggregate_downstream.py"