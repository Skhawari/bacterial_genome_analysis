# rules/quast.smk

# Assembly quality control rules using QUAST and BUSCO
# Provides comprehensive assembly metrics and completeness assessment

rule busco_completeness:
    input:
        contigs = "results/assembly/{sample}/contigs.fasta"
    output:
        summary = "results/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.txt",
        full_table = "results/busco/{sample}/run_bacteria_odb10/full_table.tsv"
    log:
        "logs/busco/{sample}.log"
    threads: 4
    resources:
        mem_mb = 8000
    conda:
        "../envs/quast.yaml"
    shell:
        """
        busco -i {input.contigs} \
              -o {wildcards.sample} \
              --out_path results/busco/ \
              --lineage_dataset bacteria_odb10 \
              --download_path busco_downloads/ \
              --mode genome \
              --cpu {threads} \
              --force \
              > {log} 2>&1
        """

rule quast_assembly_qc:
    input:
        contigs = "results/assembly/{sample}/contigs.fasta"
    output:
        report = "results/assembly_qc/{sample}/report.html",
        stats = "results/assembly_qc/{sample}/report.txt"
    log:
        "logs/quast/{sample}.log"
    threads: 4
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast.py {input.contigs} \
                 -o results/assembly_qc/{wildcards.sample} \
                 --threads {threads} \
                 --min-contig 500 \
                 --est-ref-size 5000000 \
                 > {log} 2>&1
        """

rule multiqc_assembly:
    input:
        quast_reports = expand("results/assembly_qc/{sample}/report.txt", sample=ASSEMBLY_SAMPLES),
        busco_summaries = expand("results/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.txt", sample=ASSEMBLY_SAMPLES)
    output:
        report = "results/assembly_multiqc/assembly_multiqc_report.html"
    log:
        "logs/multiqc_assembly.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc results/assembly_qc/ results/busco/ \
                -o results/assembly_multiqc/ \
                -n assembly_multiqc_report \
                --title "Assembly Quality Assessment Report" \
                --comment "Comprehensive assembly quality metrics from QUAST and BUSCO" \
                --force \
                > {log} 2>&1
        """

rule assembly_summary:
    input:
        quast_reports = expand("results/assembly_qc/{sample}/report.txt", sample=ASSEMBLY_SAMPLES),
        busco_summaries = expand("results/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.txt", sample=ASSEMBLY_SAMPLES)
    output:
        summary = "results/assembly/assembly_summary.txt"
    log:
        "logs/assembly_summary.log"
    shell:
        """
        echo "Assembly Quality Assessment Summary" > {output.summary}
        echo "=================================" >> {output.summary}
        echo "Global configuration mode: {config[assembly][mode]}" >> {output.summary}
        echo "Total samples processed: {ASSEMBLY_SAMPLES}" >> {output.summary}
        echo "Samples with long reads available: {LONG_READ_AVAILABLE}" >> {output.summary}
        echo "Samples blacklisted from long-read analysis: {LONG_READ_BLACKLISTED}" >> {output.summary}
        echo "" >> {output.summary}
        
        echo "QUAST Results Summary:" >> {output.summary}
        echo "---------------------" >> {output.summary}
        for quast_file in {input.quast_reports}; do
            sample=$(basename $(dirname $quast_file))
            echo "Sample: $sample" >> {output.summary}
            grep -E "(# contigs|Total length|N50|Largest contig)" $quast_file >> {output.summary}
            echo "" >> {output.summary}
        done
        
        echo "BUSCO Completeness Summary:" >> {output.summary}
        echo "--------------------------" >> {output.summary}
        for busco_file in {input.busco_summaries}; do
            sample=$(basename $busco_file | cut -d'.' -f4)
            echo "Sample: $sample" >> {output.summary}
            grep -E "(Complete BUSCOs|Missing BUSCOs)" $busco_file >> {output.summary}
            echo "" >> {output.summary}
        done
        """
