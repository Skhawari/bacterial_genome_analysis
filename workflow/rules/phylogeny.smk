# rules/phylogeny.smk

# Phylogenetic analysis rules using IQ-TREE
# Constructs maximum likelihood phylogenetic trees from core genome alignments

rule iqtree_phylogeny:
    input:
        alignment = "results/comparative/core_alignment.phy",
        alignment_info = "results/comparative/alignment_info.txt"
    output:
        # IQ-TREE main outputs
        tree = "results/phylogeny/core_genome.treefile",
        log = "results/phylogeny/core_genome.log",
        iqtree = "results/phylogeny/core_genome.iqtree",
        # Model selection and bootstrap outputs
        model_file = "results/phylogeny/core_genome.model.gz",
        bootstrap_trees = "results/phylogeny/core_genome.boottrees",
        consensus_tree = "results/phylogeny/core_genome.contree",
        # Additional outputs
        splits = "results/phylogeny/core_genome.splits.nex",
        checkpoint = "results/phylogeny/core_genome.ckp.gz"
    params:
        prefix = "results/phylogeny/core_genome",
        model_selection = config["phylogeny"].get("model_selection", "MFP"),
        bootstrap_replicates = config["phylogeny"].get("bootstrap_replicates", 1000),
        extra = config["phylogeny"].get("extra", ""),
        outgroup = config.get("outgroup_name", ""),
        threads_iqtree = lambda wildcards, threads: min(threads, 16)  # IQ-TREE works best with ≤16 threads
    log:
        "logs/phylogeny/iqtree_phylogeny.log"
    threads: workflow.cores if workflow.cores <= 16 else 16
    resources:
        mem_mb = 32000
    conda:
        "../envs/phylogeny.yaml"
    shell:
        """
        # Copy alignment to working directory with correct name
        cp {input.alignment} {params.prefix}.phy
        
        echo "Starting IQ-TREE phylogenetic analysis..." | tee -a {log}
        echo "Alignment: {input.alignment}" | tee -a {log}
        echo "Model selection: {params.model_selection}" | tee -a {log}
        echo "Bootstrap replicates: {params.bootstrap_replicates}" | tee -a {log}
        echo "Threads: {params.threads_iqtree}" | tee -a {log}
        echo "Memory: {resources.mem_mb}MB"  | tee -a {log}
        
        # Build IQ-TREE command
        iqtree_cmd="iqtree -s {params.prefix}.phy -m {params.model_selection} -bb {params.bootstrap_replicates} -nt {params.threads_iqtree}"
        
        # Add outgroup if specified
        if [ -n "{params.outgroup}" ]; then
            echo "Using outgroup: {params.outgroup}" | tee -a {log}
            iqtree_cmd="$iqtree_cmd -o {params.outgroup}"
        fi
        
        # Add extra parameters if specified
        if [ -n "{params.extra}" ]; then
            echo "Extra parameters: {params.extra}" | tee -a {log}
            iqtree_cmd="$iqtree_cmd {params.extra}"
        fi
        
        # Run IQ-TREE
        echo "Running command: $iqtree_cmd" | tee -a {log}
        eval $iqtree_cmd 2>&1 | tee -a {log}
        
        # Move output files to expected locations (IQ-TREE adds .phy to all output names)
        mv {params.prefix}.phy.treefile {output.tree}
        mv {params.prefix}.phy.log {output.log}
        mv {params.prefix}.phy.iqtree {output.iqtree}
        mv {params.prefix}.phy.model.gz {output.model_file}
        mv {params.prefix}.phy.contree {output.consensus_tree}
        mv {params.prefix}.phy.splits.nex {output.splits}
        mv {params.prefix}.phy.ckp.gz {output.checkpoint}
        
        # Handle optional bootstrap trees file
        if [ -f {params.prefix}.phy.boottrees ]; then
            mv {params.prefix}.phy.boottrees {output.bootstrap_trees}
        else
            touch {output.bootstrap_trees}  # Create empty file if not generated
        fi
        
        # Check if main outputs were created
        if [ ! -f {output.tree} ]; then
            echo "Error: Tree file not created" | tee -a {log}
            exit 1
        fi
        
        if [ ! -f {output.iqtree} ]; then
            echo "Error: IQ-TREE report not created" | tee -a {log}
            exit 1
        fi
        
        echo "IQ-TREE phylogenetic analysis completed successfully!" | tee -a {log}
        echo "Main tree file: {output.tree}" | tee -a {log}
        echo "Bootstrap consensus tree: {output.consensus_tree}" | tee -a {log}
        """

rule phylogeny_visualization:
    input:
        tree = "results/phylogeny/core_genome.treefile",
        consensus_tree = "results/phylogeny/core_genome.contree",
        iqtree_log = "results/phylogeny/core_genome.iqtree",
        alignment_info = "results/comparative/alignment_info.txt"
    output:
        tree_plot = "results/phylogeny/phylogenetic_tree.png",
        tree_plot_pdf = "results/phylogeny/phylogenetic_tree.pdf",
        consensus_plot = "results/phylogeny/consensus_tree.png",
        html_report = "results/phylogeny/phylogeny_report.html",
        tree_stats = "results/phylogeny/tree_statistics.txt"
    log:
        "logs/phylogeny/phylogeny_visualization.log"
    threads: 8
    conda:
        "../envs/phylogeny.yaml"
    script:
        "../scripts/phylogeny_visualization.py"

rule phylogeny_summary:
    input:
        tree = "results/phylogeny/core_genome.treefile",
        consensus_tree = "results/phylogeny/core_genome.contree",
        iqtree_log = "results/phylogeny/core_genome.iqtree",
        html_report = "results/phylogeny/phylogeny_report.html",
        tree_plot = "results/phylogeny/phylogenetic_tree.png"
    output:
        summary_report = "results/phylogeny/phylogeny_summary.txt"
    log:
        "logs/phylogeny/phylogeny_summary.log"
    shell:
        """
        echo "Phylogenetic Analysis Summary" > {output.summary_report}
        echo "=============================" >> {output.summary_report}
        echo "Generated on: $(date)" >> {output.summary_report}
        echo "" >> {output.summary_report}
        
        echo "Files generated:" >> {output.summary_report}
        echo "- Maximum likelihood tree: {input.tree}" >> {output.summary_report}
        echo "- Bootstrap consensus tree: {input.consensus_tree}" >> {output.summary_report}
        echo "- Tree visualization: {input.tree_plot}" >> {output.summary_report}
        echo "- Detailed report: {input.html_report}" >> {output.summary_report}
        echo "" >> {output.summary_report}
        
        echo "Tree files are in Newick format and can be opened with:" >> {output.summary_report}
        echo "- FigTree (https://github.com/rambaut/figtree/)" >> {output.summary_report}
        echo "- iTOL (https://itol.embl.de/)" >> {output.summary_report}
        echo "- Dendroscope" >> {output.summary_report}
        echo "- R packages: ape, ggtree, phytools" >> {output.summary_report}
        echo "" >> {output.summary_report}
        
        # Extract model information from IQ-TREE log
        if [ -f {input.iqtree_log} ]; then
            echo "Best-fit model information:" >> {output.summary_report}
            grep -A 5 "Best-fit model:" {input.iqtree_log} >> {output.summary_report} || echo "Model info not found" >> {output.summary_report}
            echo "" >> {output.summary_report}
        fi
        
        echo "Phylogenetic analysis completed successfully!" >> {output.summary_report}
        echo "Summary written to: {output.summary_report}" | tee -a {log}
        """
