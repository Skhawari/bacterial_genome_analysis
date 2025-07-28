rule aggregate_final_report:
    input:
        qc_report        = "results/multiqc/multiqc_report.html",
        assembly_qa      = "results/assembly_multiqc/assembly_multiqc_report.html",
        comp_report      = "results/comparative/comparative_genomics_report.html",
        phylo_report     = "results/phylogeny/phylogeny_report.html",
        downstream_sum   = "results/downstream/summary.xlsx" if config["pipeline"].get("run_downstream", False) else "config/config.yaml"
    output:
        html = "reports/final_report.html",
        qc_copy = "reports/multiqc_report.html",
        assembly_copy = "reports/assembly_multiqc_report.html",
        comp_copy = "reports/comparative_genomics_report.html",
        phylo_copy = "reports/phylogeny_report.html"
    params:
        title = config.get("reporting", {}).get("title", "Bacterial Genome Analysis Pipeline Report"),
        downstream_enabled = config["pipeline"].get("run_downstream", False)
    log:
        "logs/final_report.log"
    conda:
        "../envs/reporting.yaml"
    run:
        # Copy HTML files to reports directory (Snakemake creates the output directory)
        shell("cp {input.qc_report} {output.qc_copy}")
        shell("cp {input.assembly_qa} {output.assembly_copy}")
        shell("cp {input.comp_report} {output.comp_copy}")
        shell("cp {input.phylo_report} {output.phylo_copy}")
        
        # Run the Python script to generate the final report
        import pandas as pd
        import jinja2
        from pathlib import Path

        # Eingabedateien aus snakemake
        downstream_xlsx = input.downstream_sum
        title           = params.title
        enable_downstream = params.downstream_enabled

        # Versuche Downstream-Daten zu laden
        ds_dict = {}
        if enable_downstream and downstream_xlsx.endswith(".xlsx") and Path(downstream_xlsx).exists():
            try:
                ds_dict = pd.read_excel(downstream_xlsx, sheet_name=None)
            except Exception as e:
                print(f"⚠️ Warning: Could not read Excel {downstream_xlsx}: {e}")
                ds_dict = {}

        # Jinja2-Template vorbereiten
        env = jinja2.Environment(loader=jinja2.FileSystemLoader("templates"))
        tmpl = env.get_template("final_report.html.j2")

        # HTML-Rendering
        html_content = tmpl.render(
            title=title,
            qc_iframe="multiqc_report.html",
            assembly_iframe="assembly_multiqc_report.html",
            comp_iframe="comparative_genomics_report.html",
            phylo_iframe="phylogeny_report.html",
            downstream_tables={
                name: df.to_html(classes="table table-striped table-bordered", index=False, escape=True)
                for name, df in ds_dict.items()
            }
        )

        # HTML schreiben
        with open(output.html, "w") as fh:
            fh.write(html_content)