rule aggregate_final_report:
    input:
        # HTML Reports
        qc_report = "results/multiqc/multiqc_report.html",
        assembly_qa = "results/assembly_multiqc/assembly_multiqc_report.html",
        comp_report = "results/comparative/comparative_genomics_report.html",
        phylo_report = "results/phylogeny/phylogeny_report.html",
        tree_img = "results/phylogeny/phylogenetic_tree.png",
        prokka_report = "results/visualization/prokka_annotation_report_interactive.html",
        # Make downstream summary
        downstream_sum = "results/downstream/summary.xlsx"
    output:
        html = "reports/final_report.html",
        qc_copy = "reports/multiqc_report.html",
        assembly_copy = "reports/assembly_multiqc_report.html",
        comp_copy = "reports/comparative_genomics_report.html",
        phylo_copy = "reports/phylogeny_report.html",
        tree_img_copy = "reports/phylogenetic_tree.png",
        prokka_copy = "reports/prokka_report.html",
        downstream_copy = "reports/downstream_sum.xlsx"
    params:
        title = config.get("reporting", {}).get("title", "Bacterial Genome Analysis Pipeline Report"),
        downstream_enabled = any([
            config["downstream"].get("run_mlst", False),
            config["downstream"].get("run_amr", False),
            config["downstream"].get("run_virulence", False),
            config["downstream"].get("run_plasmid", False)
        ])
    log:
        "logs/final_report/final_report.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 30
    conda:
        "../envs/reporting.yaml"
    script:
        "../scripts/render_report.py"