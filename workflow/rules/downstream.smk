# rules/downstream.smk

# Conditional flags from config
run_mlst      = config["downstream"].get("run_mlst", False)
run_amr       = config["downstream"].get("run_amr", False)
run_virulence = config["downstream"].get("run_virulence", False)
run_plasmid   = config["downstream"].get("run_plasmid", False)

def get_contigs_for_downstream(wildcards):
    """Get contigs file for downstream analysis - either assembly output or outgroup FASTA"""
    outgroup_name = config.get("outgroup_name", "outgroup")
    if wildcards.sample == outgroup_name:
        return config["outgroup_fasta"]
    else:
        return f"results/assembly/{wildcards.sample}/contigs.fasta"

# MLST Typing
rule mlst_typing:
    input:
        contigs = get_contigs_for_downstream
    output:
        "results/downstream/mlst/{sample}.txt"
    params:
        scheme = config["mlst"]["scheme"],
        skip = not run_mlst
    log:
        "logs/mlst/{sample}.log"
    conda: 
        "../envs/downstream.yaml"
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        """
        if ["{params.skip}" = "True"]; then
            echo "Skipping MLST for {wildcards.sample}" >> {log}
            touch {output}
        else
            mlst --scheme {params.scheme} --threads {threads} {input.contigs} > {output} 2> {log}
        fi
        """

# AMR Screen
rule amr_screen:
    input:
        contigs = get_contigs_for_downstream
    output:
        "results/downstream/amr/{sample}.tsv"
    params:
        extra = config["amr"]["extra"],
        skip = not run_amr
    log:
        "logs/amr/{sample}.log"
    conda: 
        "../envs/downstream.yaml"
    threads: 16
    resources:
        mem_mb = 8000
    shell:
        """
        if ["{params.skip}" = "True"]; then
            echo "Skipping AMR for {wildcards.sample}" >> {log}
            touch {output}
        else
            abricate --db card --threads {threads} {params.extra} {input.contigs} > {output} 2> {log}
        fi
        """

# Virulence Screen
rule virulence_screen:
    input:
        contigs = get_contigs_for_downstream
    output:
        "results/downstream/virulence/{sample}.tsv"
    params:
        extra = config["virulence"]["extra"],
        skip = not run_virulence
    log:
        "logs/virulence/{sample}.log"
    conda: 
        "../envs/downstream.yaml"
    threads: 16
    resources:
        mem_mb = 8000
    shell:
        """
        if ["{params.skip}" = "True"]; then
            echo "Skipping virulence for {wildcards.sample}" >> {log}
            touch {output}
        else
            abricate --db vfdb --threads {threads} {params.extra} {input.contigs} > {output} 2> {log}
        fi
        """

# PLSDB Sketch
rule mash_sketch_plasmid_db:
    input:
        fasta = config["plasmid_db"]["fasta"]
    output:
        "databases/plsdb/plasmids.msh"
    log:
        "logs/plasmid/mash_sketch_db.log"
    conda: 
        "../envs/downstream.yaml"
    threads: 16
    resources:
        mem_mb = 8000
    shell:
        """
        mash sketch -p {threads} -o databases/plsdb/plasmids {input.fasta} > {log} 2>&1
        """

# PLSDB BLAST DB
rule makeblastdb_plasmid:
    input:
        fasta = config["plasmid_db"]["fasta"]
    output:
        nsq= "databases/plsdb/plasmids.nsq",
        nin= "databases/plsdb/plasmids.nin",
        nhr= "databases/plsdb/plasmids.nhr"
    log:
        "logs/plasmid/makeblastdb.log"
    conda: 
        "../envs/downstream.yaml"
    threads: 1
    shell:
        """
        makeblastdb -in {input.fasta} -dbtype nucl -out databases/plsdb/plasmids > {log} 2>&1
        """

# Plasmid Mash Screen
rule mash_screen_plasmids:
    input:
        contigs = get_contigs_for_downstream,
        db = "databases/plsdb/plasmids.msh"
    output:
        "results/downstream/plasmid/{sample}/mash_hits.txt"
    params:
        thresh = config["plasmid_db"]["mash_threshold"],
        skip = not run_plasmid
    log:
        "logs/plasmid/{sample}.mash.log"
    conda: 
        "../envs/downstream.yaml"
    threads: 16
    resources:
        mem_mb = 8000
    shell:
        """
        if ["{params.skip}" = "True"]; then
            echo "Skipping plasmid screen" >> {log}
            touch {output}
        else
            mash screen -p {threads} -i {params.thresh} {input.db} {input.contigs} > {output} 2> {log}
        fi
        """

# Plasmid BLASTn
rule plasmid_blastn:
    input:
        contigs    = get_contigs_for_downstream,
        hits       = "results/downstream/plasmid/{sample}/mash_hits.txt",
        nsq        = "databases/plsdb/plasmids.nsq"
    output:
        "results/downstream/plasmid/{sample}/blast_plasmid_hits.tsv"
    params:
        blast_db = "databases/plsdb/plasmids",
        evalue = config["plasmid_db"]["blast_evalue"],
        pid    = config["plasmid_db"]["blast_pid"],
        skip = not run_plasmid
    log:
        "logs/plasmid/{sample}.blast.log"
    conda: 
        "../envs/downstream.yaml"
    threads: 32
    resources:
        mem_mb = 16000
    script:
        "../scripts/plasmid_blast.py"

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
    input:
        entries = define_agg_targets()
    output:
        "results/downstream/summary.xlsx"
    log:
        "logs/downstream/aggregate.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: 
        "../envs/downstream.yaml"
    script:
        "../scripts/aggregate_downstream.py"