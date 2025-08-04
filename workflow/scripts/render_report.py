#!/usr/bin/env python3
import pandas as pd
import jinja2
from pathlib import Path
import shutil
from datetime import datetime

def main():
    log_file = snakemake.log[0]
    
    try:
        with open(log_file, 'w') as f:
            f.write("Starting final report generation...\n")
        
        downstream_tables = {}

        # Copy HTML reports
        shutil.copy2(snakemake.input.qc_report, snakemake.output.qc_copy)
        shutil.copy2(snakemake.input.assembly_qa, snakemake.output.assembly_copy)
        shutil.copy2(snakemake.input.comp_report, snakemake.output.comp_copy)
        shutil.copy2(snakemake.input.phylo_report, snakemake.output.phylo_copy)
        shutil.copy2(snakemake.input.tree_img,       snakemake.output.tree_img_copy)
        shutil.copy2(snakemake.input.prokka_report,  snakemake.output.prokka_copy)
        shutil.copy2(snakemake.input.downstream_sum, snakemake.output.downstream_copy)
        
        # Check if downstream is enabled AND file exists
        if snakemake.params.downstream_enabled:
            if hasattr(snakemake.input, 'downstream_sum') and snakemake.input.downstream_sum:
                downstream_file = snakemake.input.downstream_sum
                if isinstance(downstream_file, list) and len(downstream_file) > 0:
                    downstream_file = downstream_file[0]
                
                if downstream_file and Path(downstream_file).exists():
                    try:
                        all_sheets = pd.read_excel(downstream_file, sheet_name=None)
                        
                        with open(log_file, 'a') as f:
                            f.write(f"Loaded downstream Excel with {len(all_sheets)} sheets: {list(all_sheets.keys())}\n")
                        
                        for sheet_name, df in all_sheets.items():
                            if not df.empty:
                                table_html = df.to_html(
                                    classes="table table-striped table-bordered table-hover",
                                    index=False,
                                    escape=False
                                )
                                downstream_tables[f"{sheet_name} Results"] = table_html
                            
                    except Exception as e:
                        with open(log_file, 'a') as f:
                            f.write(f"Error loading downstream Excel: {e}\n")
                else:
                    with open(log_file, 'a') as f:
                        f.write("Downstream enabled but file not found\n")
            else:
                with open(log_file, 'a') as f:
                    f.write("Downstream enabled but no input file specified\n")
        else:
            with open(log_file, 'a') as f:
                f.write("Downstream analysis disabled\n")
        
        # === 4. Render template ===
        env = jinja2.Environment(loader=jinja2.FileSystemLoader("workflow/templates"))
        template = env.get_template("final_report.html.j2")
        
        html_content = template.render(
            title=snakemake.params.title,
            qc_iframe="multiqc_report.html",
            assembly_iframe="assembly_multiqc_report.html",
            comp_iframe="comparative_genomics_report.html",
            phylo_iframe="phylogeny_report.html",
            tree_image="phylogenetic_tree.png",
            prokka_iframe="prokka_report.html",
            
            # Downstream - Handle both enabled/disabled cases
            downstream_tables=downstream_tables,
            has_downstream_results=len(downstream_tables) > 0,
            downstream_enabled=snakemake.params.downstream_enabled,
            now=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )
        
        # === 5. Write final HTML ===
        with open(snakemake.output.html, "w") as f:
            f.write(html_content)
        
        with open(log_file, 'a') as f:
            f.write(f"Final report generated successfully!\n")
            f.write(f"Downstream enabled: {snakemake.params.downstream_enabled}\n")
            f.write(f"Downstream sections: {len(downstream_tables)}\n")
            
    except Exception as e:
        with open(log_file, 'a') as f:
            f.write(f"ERROR: {e}\n")
        raise

if __name__ == "__main__":
    main()