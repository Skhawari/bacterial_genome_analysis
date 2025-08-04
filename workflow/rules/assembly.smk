rule spades_hybrid:
    input:
        fq1  = "results/trimmed/{sample}_R1.fastq.gz",
        fq2  = "results/trimmed/{sample}_R2.fastq.gz"
    output:
        contigs       = "results/assembly/{sample}/contigs.fasta",
        scaffolds     = "results/assembly/{sample}/scaffolds.fasta",
        assembly_graph= "results/assembly/{sample}/assembly_graph.fastg",
        log_file      = "results/assembly/{sample}/spades.log",
        # tmpdir als temporäres Verzeichnis
        tmpdir        = temp(directory("results/assembly/{sample}/tmp"))
    params:
        outdir    = "results/assembly/{sample}",
        mode      = lambda wc: get_assembly_mode_for_sample(wc.sample),
        long_reads= lambda wc: (
            f"results/filtered_long_reads/{wc.sample}.fastq.gz"
            if get_assembly_mode_for_sample(wc.sample) == "hybrid"
            else ""
        )
    log:
        "logs/assembly/{sample}.log"
    threads: 16
    resources:
        mem_mb = 32000
    conda:
        "../envs/assembly.yaml"
    shell:
        r"""

        if [ "{params.mode}" = "short" ]; then
            spades.py \
              -1 {input.fq1} -2 {input.fq2} \
              -o {params.outdir} \
              --tmp-dir {output.tmpdir} \
              --threads {threads} \
              --memory $(({resources.mem_mb} / 1024)) \
              {config[assembly][spades][extra]} \
              2>&1 | tee -a {log}
        else
            spades.py \
              -1 {input.fq1} -2 {input.fq2} \
              --nanopore {params.long_reads} \
              -o {params.outdir} \
              --tmp-dir {output.tmpdir} \
              --threads {threads} \
              --memory $(({resources.mem_mb} / 1024)) \
              {config[assembly][spades][extra]} \
              2>&1 | tee -a {log}
        fi
        """