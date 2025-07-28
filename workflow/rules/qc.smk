# rules/qc.smk

# Short-read QC rules
rule fastqc_raw:
    input:
        lambda wc: samples.at[wc.sample, f"fq{wc.idx}"]
    output:
        html = "results/fastqc_raw/{sample}_{idx}_fastqc.html",
        zip = "results/fastqc_raw/{sample}_{idx}_fastqc.zip"
    log:
        "logs/fastqc_raw/{sample}_{idx}.log"
    threads: 4
    resources:
        mem_mb = 1000
    wrapper:
        "v6.2.0/bio/fastqc"

rule fastp:
    input:
        sample = lambda wc: [samples.at[wc.sample, "fq1"], samples.at[wc.sample, "fq2"]]
    output:
        trimmed = ["results/trimmed/{sample}_R1.fastq.gz", "results/trimmed/{sample}_R2.fastq.gz"],
        html = "results/qc/fastp/{sample}.html",
        json = "results/qc/fastp/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 8
    wrapper:
        "v6.2.0/bio/fastp"

rule fastqc_trimmed:
    input:
        "results/trimmed/{sample}_R{idx}.fastq.gz"
    output:
        html = "results/fastqc_trimmed/{sample}_R{idx}_fastqc.html",
        zip = "results/fastqc_trimmed/{sample}_R{idx}_fastqc.zip"
    log:
        "logs/fastqc_trimmed/{sample}_{idx}.log"
    threads: 4
    resources:
        mem_mb = 1000
    wrapper:
        "v6.2.0/bio/fastqc"

# Long-read QC rules
rule nanostat_raw:
    input:
        fastq = lambda wc: long_read_samples.at[wc.sample, "long_read"]
    output:
        stats = "results/nanostat_raw/{sample}.txt"
    log:
        "logs/nanostat_raw/{sample}.log"
    threads: 1
    conda:
        "../envs/qc.yaml"
    shell:
        """
        NanoStat --fastq {input.fastq} > {output.stats} 2> {log}
        """

rule nanostat_filtered:
    input:
        fastq = "results/filtered_long_reads/{sample}.fastq.gz"
    output:
        stats = "results/nanostat_filtered/{sample}.txt"
    log:
        "logs/nanostat_filtered/{sample}.log"
    threads: 1
    conda:
        "../envs/qc.yaml"
    shell:
        """
        NanoStat --fastq {input.fastq} > {output.stats} 2> {log}
        """

# Optional: Long-read filtering with Filtlong (for bacterial genomes, typically filter reads >1kb and keep best quality)
rule filtlong_filter:
    input:
        fastq = lambda wc: long_read_samples.at[wc.sample, "long_read"]
    output:
        filtered = "results/filtered_long_reads/{sample}.fastq.gz"
    log:
        "logs/filtlong/{sample}.log"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        filtlong --min_length {config[filtlong][min_length]} \
                 --keep_percent {config[filtlong][keep_percent]} \
                 --target_bases {config[filtlong][target_bases]} \
                 {input.fastq} | gzip > {output.filtered} 2> {log}
        """

rule multiqc:
    input:
        # Short-read QC files (always present for all samples)
        expand("results/fastqc_raw/{sample}_{idx}_fastqc.zip", sample=ALL_SAMPLES, idx=["1", "2"]),
        expand("results/fastqc_trimmed/{sample}_R{idx}_fastqc.zip", sample=ALL_SAMPLES, idx=["1", "2"]),
        expand("results/qc/fastp/{sample}.json", sample=ALL_SAMPLES),
        # Long-read QC files (only for samples with long reads available)
        expand("results/nanostat_raw/{sample}.txt", sample=LONG_READ_AVAILABLE),
        expand("results/nanostat_filtered/{sample}.txt", sample=LONG_READ_AVAILABLE)
    output:
        html = "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc --dirs {input} -o results/multiqc --force > {log} 2>&1
        """