# rules/visualization.smk

# Visualization rules for Prokka annotation results
# Creates comprehensive reports and plots for annotation analysis

rule prokka_visualization:
    input:
        gff_files = expand("results/annotation/{sample}/{sample}.gff", sample=ASSEMBLY_SAMPLES),
        txt_files = expand("results/annotation/{sample}/{sample}.txt", sample=ASSEMBLY_SAMPLES),
        tsv_files = expand("results/annotation/{sample}/{sample}.tsv", sample=ASSEMBLY_SAMPLES)
    output:
        html_report = "results/visualization/prokka_annotation_report.html",
        pdf_report = "results/visualization/prokka_annotation_report.pdf",
        summary_plot = "results/visualization/annotation_summary.png",
        functional_plot = "results/visualization/functional_categories.png"
    log:
        "logs/prokka_visualization.log"
    threads: 4
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/prokka_visualization.py"

rule enhanced_data_visualization:
    input:
        # Assembly quality reports
        quast_reports = expand("results/assembly_qc/{sample}/report.txt", sample=ASSEMBLY_SAMPLES),
        # BUSCO completeness reports
        busco_reports = expand("results/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.txt", sample=ASSEMBLY_SAMPLES),
        # GFF files for gene counts (this ensures Prokka has completed)
        gff_files = expand("results/annotation/{sample}/{sample}.gff", sample=ASSEMBLY_SAMPLES),
        # Assembly summary for additional plots
        assembly_summary = "results/assembly/assembly_summary.txt",
        # Configuration file
        config_file = "config/config.yaml"
    output:
        # Comprehensive report combining all statistical analysis
        comprehensive_report = "results/enhanced_visualization/comprehensive_analysis_report.pdf"
    log:
        "logs/enhanced_data_visualization.log"
    threads: 4
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/enhanced_data_visualization.py"
