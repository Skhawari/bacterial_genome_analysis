# rules/comparative.smk

# Comparative genomics rules for core genome analysis and phylogeny
# Uses Roary for pan-genome analysis and prepares for IQ-TREE phylogeny

rule roary_pangenome:
    input:
        gff_files = expand("results/annotation/{sample}/{sample}.gff", sample=ASSEMBLY_SAMPLES),
        outgroup_gff = lambda wildcards: config["comparative"]["outgroup_gff"] if config["comparative"].get("outgroup_gff") else []
    output:
        # Core genome outputs
        core_alignment = "results/comparative/core_gene_alignment.aln",
        gene_presence_absence = "results/comparative/gene_presence_absence.csv",
        gene_presence_absence_rtab = "results/comparative/gene_presence_absence.Rtab",
        summary_statistics = "results/comparative/summary_statistics.txt",
        # Pan-genome analysis outputs
        accessory_binary_genes = "results/comparative/accessory_binary_genes.fa.newick",
        core_accessory_graph = "results/comparative/core_accessory_graph.dot",
        number_of_conserved_genes = "results/comparative/number_of_conserved_genes.Rtab",
        number_of_genes_in_pan_genome = "results/comparative/number_of_genes_in_pan_genome.Rtab",
        number_of_new_genes = "results/comparative/number_of_new_genes.Rtab",
        number_of_unique_genes = "results/comparative/number_of_unique_genes.Rtab",
        pan_genome_reference = "results/comparative/pan_genome_reference.fa",
        clustered_proteins = "results/comparative/clustered_proteins"
    params:
        output_dir = "results/comparative",
        identity_threshold = config["comparative"].get("identity_threshold", 95),
        core_definition = config["comparative"].get("core_definition", 99),
        extra = config["comparative"].get("extra", "-e -n -v"),
        outgroup_name = config["comparative"].get("outgroup_name", "outgroup")
    log:
        "logs/roary_pangenome.log"
    threads: workflow.cores if workflow.cores <= 16 else 16  # Use available cores, max 16
    conda:
        "../envs/comparative.yaml"
    shell:
        """
        # Create temporary directory for input files
        mkdir -p temp_gff_roary
        
        # Copy sample GFF files with proper naming
        for gff in {input.gff_files}; do
            sample=$(basename $gff .gff)
            cp $gff temp_gff_roary/$sample.gff
        done
        
        # Copy outgroup if provided (check if outgroup_gff is not empty)
        if [ -n "{input.outgroup_gff}" ] && [ "{input.outgroup_gff}" != "[]" ]; then
            echo "Including outgroup: {input.outgroup_gff} as {params.outgroup_name}" | tee -a {log}
            cp {input.outgroup_gff} temp_gff_roary/{params.outgroup_name}.gff
        else
            echo "No outgroup provided - running analysis on samples only" | tee -a {log}
        fi
        
        # Count input files
        n_files=$(ls temp_gff_roary/*.gff | wc -l)
        echo "Running Roary on $n_files genomes" | tee -a {log}
        
        # Run Roary pan-genome analysis
        roary -f {params.output_dir} \
              -p {threads} \
              -i {params.identity_threshold} \
              -cd {params.core_definition} \
              {params.extra} \
              temp_gff_roary/*.gff 2>&1 | tee -a {log}
        
        # Clean up temporary directory
        rm -rf temp_gff_roary
        
        # Roary creates output in a directory with random suffix, so we need to find and move files
        roary_output_dir=$(find results/ -name "comparative_*" -type d | head -1)
        
        if [ -n "$roary_output_dir" ] && [ -d "$roary_output_dir" ]; then
            echo "Found Roary output directory: $roary_output_dir" | tee -a {log}
            
            # Move all files to the expected directory
            mkdir -p {params.output_dir}
            mv "$roary_output_dir"/* {params.output_dir}/
            rmdir "$roary_output_dir"
            
            echo "Successfully moved Roary outputs to {params.output_dir}" | tee -a {log}
        else
            echo "Error: Could not find Roary output directory" | tee -a {log}
            exit 1
        fi
        
        # Verify that required files exist
        for file in {output.gene_presence_absence} {output.core_alignment}; do
            if [ ! -f "$file" ]; then
                echo "Error: Required output file $file not found" | tee -a {log}
                exit 1
            fi
        done
        
        echo "Roary pan-genome analysis completed successfully!" | tee -a {log}
        """

rule roary_visualization:
    input:
        gene_presence_absence = "results/comparative/gene_presence_absence.csv",
        gene_presence_absence_rtab = "results/comparative/gene_presence_absence.Rtab",
        number_of_conserved_genes = "results/comparative/number_of_conserved_genes.Rtab",
        number_of_genes_in_pan_genome = "results/comparative/number_of_genes_in_pan_genome.Rtab"
    output:
        pangenome_plot = "results/comparative/pangenome_analysis.png",
        core_accessory_plot = "results/comparative/core_accessory_distribution.png",
        gene_frequency_plot = "results/comparative/gene_frequency_distribution.png",
        html_report = "results/comparative/comparative_genomics_report.html"
    log:
        "logs/roary_visualization.log"
    threads: 2
    conda:
        "../envs/comparative.yaml"
    script:
        "../scripts/comparative_visualization.py"

rule prepare_phylogeny_input:
    input:
        core_alignment = "results/comparative/core_gene_alignment.aln"
    output:
        phylip_alignment = "results/comparative/core_alignment.phy",
        alignment_info = "results/comparative/alignment_info.txt"
    log:
        "logs/prepare_phylogeny.log"
    conda:
        "../envs/comparative.yaml"
    script:
        "../scripts/prepare_phylogeny.py"
