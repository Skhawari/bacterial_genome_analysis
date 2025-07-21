# rules/reporting.smk

configfile: "config/config.yaml"
import pandas as pd

# Samples ohne Blacklist
ALL_SAMPLES   = pd.read_csv(config["samples"], sep="\t").sample_id.tolist()
DROP_ALL      = set(open(config["complete_blacklist"]).read().splitlines())
ANNOT_SAMPLES = [s for s in ALL_SAMPLES if s not in DROP_ALL]

rule aggregate_final_report:
    """
    Erstellt einen Master-Report, der:
      - QC MultiQC
      - Assembly-QA MultiQC
      - Comparative Genomics HTML
      - Phylogenie HTML
      - Downstream Excel (tabellarisch eingebettet)
    zusammenfasst.
    """
    conda: "envs/reporting.yaml"
    input:
        qc_report        = "qc/multiqc/qc.html",
        assembly_qa      = "qc/multiqc/assembly_qa.html",
        comp_report      = "comparative/comparative_genomics_report.html",
        phylo_report     = "phylogeny/phylogeny_report.html",
        downstream_sum   = "downstream/summary.xlsx"
    output:
        html = "reports/final_report.html"
    params:
        title = config.get("reporting", {}).get("title", "Pipeline Summary Report")
    script:
        "scripts/render_report.py"