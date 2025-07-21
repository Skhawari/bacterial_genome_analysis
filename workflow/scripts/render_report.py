import pandas as pd
import jinja2
import sys

# Config laden
# (optional, falls du den title aus config übernehmen möchtest)

# Eingaben
qc_html        = snakemake.input.qc_report
assembly_html  = snakemake.input.assembly_qa
comp_html      = snakemake.input.comp_report
phylo_html     = snakemake.input.phylo_report
downstream_xlsx= snakemake.input.downstream_sum

# Excel in DataFrame laden
df_ds = pd.read_excel(downstream_xlsx, sheet_name=None)  # dict von DataFrames

# Template laden
env = jinja2.Environment(loader=jinja2.FileSystemLoader("templates"))
tmpl = env.get_template("final_report.html.j2")

# Rendern
html = tmpl.render(
    title=snakemake.params.title,
    qc_iframe=qc_html,
    assembly_iframe=assembly_html,
    comp_iframe=comp_html,
    phylo_iframe=phylo_html,
    downstream_tables=df_ds  # dict[str, DataFrame]
)

with open(snakemake.output.html, "w") as fo:
    fo.write(html)