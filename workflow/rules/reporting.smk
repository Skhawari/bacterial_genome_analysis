rule aggregate_final_report:
    input:
        qc_report        = "results/multiqc/multiqc_report.html",
        assembly_qa      = "results/assembly_multiqc/assembly_multiqc_report.html",
        comp_report      = "results/comparative/comparative_genomics_report.html",
        phylo_report     = "results/phylogeny/phylogeny_report.html",
        downstream_sum   = "results/downstream/summary.xlsx" if config["pipeline"].get("run_downstream", False) else "config/config.yaml"
    output:
        html = "reports/final_report.html"
    params:
        title = config.get("reporting", {}).get("title", "Bacterial Genome Analysis Pipeline Report"),
        downstream_enabled = config["pipeline"].get("run_downstream", False)
    log:
        "logs/final_report.log"
    conda:
        "../envs/reporting.yaml"
    script:
        "../scripts/render_report.py"