# rules/assembly.smk

# SPAdes assembly rules for bacterial genomes
# Supports both hybrid (short + long reads) and short-read only assemblies

rule spades_hybrid:
    input:
        fq1 = "results/trimmed/{sample}_R1.fastq.gz",
        fq2 = "results/trimmed/{sample}_R2.fastq.gz"
    output:
        contigs = "results/assembly/{sample}/contigs.fasta",
        scaffolds = "results/assembly/{sample}/scaffolds.fasta",
        assembly_graph = "results/assembly/{sample}/assembly_graph.fastg",
        log_file = "results/assembly/{sample}/spades.log"
    params:
        outdir = "results/assembly/{sample}",
        # Use the helper function to determine assembly mode
        mode = lambda wc: get_assembly_mode_for_sample(wc.sample),
        long_reads = lambda wc: (
            f"results/filtered_long_reads/{wc.sample}.fastq.gz" 
            if get_assembly_mode_for_sample(wc.sample) == "hybrid"
            else ""
        ),
    log:
        "logs/assembly/{sample}.log"
    threads: 16
    resources:
        mem_mb = 32000  
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        echo "Assembly configuration for {wildcards.sample}:" > {log}
        echo "  Global config mode: {config[assembly][mode]}" >> {log}
        echo "  SPAdes extra options: {config[assembly][spades][extra]}" >> {log}
        echo "  Sample in long_read_samples: $(if echo '{LONG_READ_AVAILABLE}' | grep -q '{wildcards.sample}'; then echo "yes"; else echo "no"; fi)" >> {log}
        echo "  Sample in long_read_blacklist: $(if echo '{LONG_READ_BLACKLISTED}' | grep -q '{wildcards.sample}'; then echo "yes"; else echo "no"; fi)" >> {log}
        echo "  Effective assembly mode: {params.mode}" >> {log}
        echo "  Long reads file: {params.long_reads}" >> {log}
        echo "" >> {log}
        
        if [ "{params.mode}" = "short" ]; then
            echo "Running SPAdes in SHORT-READ ONLY mode for {wildcards.sample}" >> {log}
            
            # Determine and log the reason for short-read only mode
            if [ "{config[assembly][mode]}" = "short" ]; then
                echo "Reason: Global config set to short-read only mode" >> {log}
            elif echo '{LONG_READ_BLACKLISTED}' | grep -q '{wildcards.sample}'; then
                echo "Reason: Sample is blacklisted from long-read analysis" >> {log}
            else
                echo "Reason: No long reads available for this sample" >> {log}
            fi
            spades.py -1 {input.fq1} -2 {input.fq2} \
                      -o {params.outdir} \
                      {config[assembly][spades][extra]} \
                      --threads {threads} --memory $(({resources.mem_mb}/1024)) \
                      >> {log} 2>&1
        else
            echo "Running SPAdes in HYBRID mode for {wildcards.sample}" >> {log}
            spades.py -1 {input.fq1} -2 {input.fq2} \
                      --nanopore {params.long_reads} \
                      -o {params.outdir} \
                      {config[assembly][spades][extra]} \
                      --threads {threads} --memory $(({resources.mem_mb}/1024)) \
                      >> {log} 2>&1
        fi
        
        echo "Assembly completed for {wildcards.sample}" >> {log}
        """
