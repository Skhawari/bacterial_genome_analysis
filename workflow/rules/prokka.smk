# rules/prokka.smk

# Genome annotation rules using Prokka
# Provides comprehensive prokaryotic genome annotation including GTF/GFF files

def get_annotation_input(wildcards):
    """Get input file for annotation - either assembly output or outgroup FASTA"""
    outgroup_name = config.get("outgroup_name", "outgroup")
    if wildcards.sample == outgroup_name:
        return config["outgroup_fasta"]
    else:
        return f"results/assembly/{wildcards.sample}/contigs.fasta"

rule prokka_annotation:
    input:
        contigs = get_annotation_input
    output:
        gff = "results/annotation/{sample}/{sample}.gff",
        gtf = "results/annotation/{sample}/{sample}.gtf", 
        faa = "results/annotation/{sample}/{sample}.faa",
        ffn = "results/annotation/{sample}/{sample}.ffn",
        fna = "results/annotation/{sample}/{sample}.fna",
        txt = "results/annotation/{sample}/{sample}.txt",
        tsv = "results/annotation/{sample}/{sample}.tsv"
    params:
        outdir = "results/annotation/{sample}"
    log:
        "logs/annotation/{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        # Remove output directory if it exists to avoid conflicts
        rm -rf {params.outdir}
        
        prokka --outdir {params.outdir} \
               --prefix {wildcards.sample} \
               --genus {config[annotation][genus]} \
               --species {config[annotation][species]} \
               --kingdom {config[annotation][kingdom]} \
               --centre "{config[annotation][centre]}" \
               --locustag {wildcards.sample} \
               --mincontiglen {config[annotation][mincontiglen]} \
               --evalue {config[annotation][evalue]} \
               --cpus {threads} \
               --force \
               --compliant \
               {input.contigs} \
               > {log} 2>&1
        
        # Convert GFF to GTF format (Prokka doesn't generate GTF by default)
        gffread {output.gff} -T -o {output.gtf} 2>> {log}
        """
