import pandas as pd
import jinja2
from pathlib import Path

# Eingabedateien aus snakemake
qc_html         = snakemake.input.qc_report
assembly_html   = snakemake.input.assembly_qa
comp_html       = snakemake.input.comp_report
phylo_html      = snakemake.input.phylo_report
downstream_xlsx = snakemake.input.downstream_sum
title           = snakemake.params.title
enable_downstream = snakemake.params.downstream_enabled

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
with open(snakemake.output.html, "w") as fh:
    fh.write(html_content)