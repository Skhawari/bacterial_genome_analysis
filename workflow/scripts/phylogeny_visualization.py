#!/usr/bin/env python3
"""
Phylogenetic Tree Visualization Script
Creates comprehensive visualizations and reports for IQ-TREE phylogenetic analysis
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from pathlib import Path
import re
from jinja2 import Template
from ete3 import Tree, TreeStyle, TextFace, CircleFace, AttrFace
import warnings
warnings.filterwarnings('ignore')

def parse_iqtree_log(iqtree_file):
    """Parse IQ-TREE output file for statistics and model information"""
    stats = {}
    
    try:
        with open(iqtree_file, 'r') as f:
            content = f.read()
            
        # Extract basic information
        if "Total number of sequences:" in content:
            match = re.search(r"Total number of sequences:\s+(\d+)", content)
            if match:
                stats['n_sequences'] = int(match.group(1))
                
        if "Total number of sites:" in content:
            match = re.search(r"Total number of sites:\s+([\d,]+)", content)
            if match:
                stats['n_sites'] = int(match.group(1).replace(',', ''))
                
        # Extract best model
        if "Best-fit model:" in content:
            match = re.search(r"Best-fit model:\s+(\S+)", content)
            if match:
                stats['best_model'] = match.group(1)
                
        # Extract log likelihood
        if "Log-likelihood of the tree:" in content:
            match = re.search(r"Log-likelihood of the tree:\s+([-\d.]+)", content)
            if match:
                stats['log_likelihood'] = float(match.group(1))
                
        # Extract AIC/BIC scores
        if "Akaike information criterion (AIC) score:" in content:
            match = re.search(r"Akaike information criterion \(AIC\) score:\s+([\d.]+)", content)
            if match:
                stats['aic'] = float(match.group(1))
                
        if "Bayesian information criterion (BIC) score:" in content:
            match = re.search(r"Bayesian information criterion \(BIC\) score:\s+([\d.]+)", content)
            if match:
                stats['bic'] = float(match.group(1))
                
        # Extract time information
        if "Total CPU time used:" in content:
            match = re.search(r"Total CPU time used:\s+([\d.]+)", content)
            if match:
                stats['cpu_time'] = float(match.group(1))
                
        if "Total wall-clock time used:" in content:
            match = re.search(r"Total wall-clock time used:\s+([\d.]+)", content)
            if match:
                stats['wall_time'] = float(match.group(1))
                
    except Exception as e:
        print(f"Error parsing IQ-TREE log: {e}")
        
    return stats

def visualize_tree_matplotlib(tree_file, output_file, title="Phylogenetic Tree"):
    """Create a simple tree visualization using matplotlib"""
    try:
        # Read tree
        tree = Tree(tree_file)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Get tree coordinates
        def get_tree_coords(node, x=0, y_positions=None, y_spacing=1):
            if y_positions is None:
                y_positions = {}
                
            if node.is_leaf():
                y_positions[node] = len(y_positions) * y_spacing
                return [(x, y_positions[node])], y_positions[node]
            else:
                child_coords = []
                child_ys = []
                for child in node.children:
                    coords, child_y = get_tree_coords(child, x + child.dist, y_positions, y_spacing)
                    child_coords.extend(coords)
                    child_ys.append(child_y)
                
                # Position internal node at midpoint of children
                y_pos = np.mean(child_ys)
                child_coords.append((x, y_pos))
                
                # Draw horizontal lines to children
                for child_y in child_ys:
                    ax.plot([x, x + min(child.dist for child in node.children)], 
                           [y_pos, child_y], 'k-', linewidth=1)
                
                return child_coords, y_pos
        
        # Get coordinates and draw tree
        coords, root_y = get_tree_coords(tree)
        
        # Plot tree structure
        for i, node in enumerate(tree.traverse("postorder")):
            if not node.is_leaf():
                continue
            # Plot leaf names
            y_pos = i
            x_pos = tree.get_distance(node)
            ax.text(x_pos + 0.001, y_pos, node.name, 
                   verticalalignment='center', fontsize=10)
        
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Evolutionary Distance')
        ax.set_ylabel('Taxa')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return True
        
    except Exception as e:
        print(f"Error creating matplotlib tree visualization: {e}")
        
        # Create a simple placeholder plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, f'Tree visualization error:\n{str(e)}\n\nPlease use external tools like FigTree\nto visualize the tree file.',
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return False

def create_tree_statistics_plot(stats, output_file):
    """Create a visualization of tree statistics"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Phylogenetic Analysis Statistics', fontsize=16, fontweight='bold')
    
    # Plot 1: Basic statistics
    if 'n_sequences' in stats and 'n_sites' in stats:
        categories = ['Sequences', 'Alignment Sites']
        values = [stats['n_sequences'], stats['n_sites']]
        colors = ['#2E86AB', '#A23B72']
        
        bars = ax1.bar(categories, values, color=colors, alpha=0.7, edgecolor='black')
        ax1.set_title('Dataset Statistics', fontweight='bold')
        ax1.set_ylabel('Count')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:,}', ha='center', va='bottom', fontweight='bold')
    else:
        ax1.text(0.5, 0.5, 'Dataset statistics\nnot available', 
                ha='center', va='center', transform=ax1.transAxes)
        ax1.set_title('Dataset Statistics', fontweight='bold')
    
    # Plot 2: Model selection scores
    if 'aic' in stats and 'bic' in stats:
        scores = ['AIC', 'BIC']
        values = [stats['aic'], stats['bic']]
        colors = ['#F18F01', '#C73E1D']
        
        bars = ax2.bar(scores, values, color=colors, alpha=0.7, edgecolor='black')
        ax2.set_title('Model Selection Scores', fontweight='bold')
        ax2.set_ylabel('Score')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:.1f}', ha='center', va='bottom', fontweight='bold')
    else:
        ax2.text(0.5, 0.5, 'Model selection scores\nnot available', 
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Model Selection Scores', fontweight='bold')
    
    # Plot 3: Best model
    if 'best_model' in stats:
        ax3.text(0.5, 0.5, f"Best-fit Model:\n{stats['best_model']}", 
                ha='center', va='center', transform=ax3.transAxes,
                fontsize=14, fontweight='bold', 
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.7))
        ax3.set_title('Model Selection Result', fontweight='bold')
        ax3.axis('off')
    else:
        ax3.text(0.5, 0.5, 'Best model\nnot available', 
                ha='center', va='center', transform=ax3.transAxes)
        ax3.set_title('Model Selection Result', fontweight='bold')
        ax3.axis('off')
    
    # Plot 4: Runtime statistics
    if 'cpu_time' in stats and 'wall_time' in stats:
        times = ['CPU Time', 'Wall Time']
        values = [stats['cpu_time'], stats['wall_time']]
        colors = ['#3D5A80', '#98C1D9']
        
        bars = ax4.bar(times, values, color=colors, alpha=0.7, edgecolor='black')
        ax4.set_title('Runtime Statistics', fontweight='bold')
        ax4.set_ylabel('Time (seconds)')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:.1f}s', ha='center', va='bottom', fontweight='bold')
    else:
        ax4.text(0.5, 0.5, 'Runtime statistics\nnot available', 
                ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('Runtime Statistics', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def create_html_report(stats, tree_files, output_file):
    """Create comprehensive HTML report for phylogenetic analysis"""
    
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Phylogenetic Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
            .header { text-align: center; color: #2E86AB; margin-bottom: 30px; }
            .section { margin: 30px 0; }
            .stats-table { width: 100%; border-collapse: collapse; margin: 20px 0; }
            .stats-table th, .stats-table td { border: 1px solid #ddd; padding: 10px; text-align: left; }
            .stats-table th { background-color: #f2f2f2; font-weight: bold; }
            .metric { display: inline-block; margin: 15px; padding: 20px; background: #f8f9fa; border-radius: 8px; }
            .metric-value { font-size: 24px; font-weight: bold; color: #2E86AB; }
            .metric-label { font-size: 14px; color: #666; margin-top: 5px; }
            .file-list { background: #f8f9fa; padding: 15px; border-radius: 8px; }
            .file-list ul { margin: 0; padding-left: 20px; }
            .model-box { background: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Phylogenetic Analysis Report</h1>
            <h2>Maximum Likelihood Tree Construction</h2>
            <p>Generated on {{ date }}</p>
        </div>
        
        <div class="section">
            <h2>Analysis Overview</h2>
            <div style="text-align: center;">
                {% if stats.n_sequences %}
                <div class="metric">
                    <div class="metric-value">{{ stats.n_sequences }}</div>
                    <div class="metric-label">Sequences</div>
                </div>
                {% endif %}
                {% if stats.n_sites %}
                <div class="metric">
                    <div class="metric-value">{{ "{:,}".format(stats.n_sites) }}</div>
                    <div class="metric-label">Alignment Sites</div>
                </div>
                {% endif %}
                {% if stats.best_model %}
                <div class="metric">
                    <div class="metric-value">{{ stats.best_model }}</div>
                    <div class="metric-label">Best Model</div>
                </div>
                {% endif %}
                {% if stats.cpu_time %}
                <div class="metric">
                    <div class="metric-value">{{ "%.1f"|format(stats.cpu_time) }}s</div>
                    <div class="metric-label">CPU Time</div>
                </div>
                {% endif %}
            </div>
        </div>
        
        {% if stats.best_model %}
        <div class="section">
            <h2>Model Selection</h2>
            <div class="model-box">
                <h3>Best-fit Model: {{ stats.best_model }}</h3>
                <p>IQ-TREE performed automatic model selection to identify the evolutionary model that best fits your data.</p>
                {% if stats.aic and stats.bic %}
                <table class="stats-table" style="width: 50%;">
                    <tr>
                        <th>Criterion</th>
                        <th>Score</th>
                    </tr>
                    <tr>
                        <td>AIC (Akaike Information Criterion)</td>
                        <td>{{ "%.2f"|format(stats.aic) }}</td>
                    </tr>
                    <tr>
                        <td>BIC (Bayesian Information Criterion)</td>
                        <td>{{ "%.2f"|format(stats.bic) }}</td>
                    </tr>
                    {% if stats.log_likelihood %}
                    <tr>
                        <td>Log-likelihood</td>
                        <td>{{ "%.2f"|format(stats.log_likelihood) }}</td>
                    </tr>
                    {% endif %}
                </table>
                {% endif %}
            </div>
        </div>
        {% endif %}
        
        <div class="section">
            <h2>Output Files</h2>
            <div class="file-list">
                <h3>Generated Files:</h3>
                <ul>
                    <li><strong>{{ tree_files.tree }}</strong> - Maximum likelihood tree (Newick format)</li>
                    <li><strong>{{ tree_files.consensus }}</strong> - Bootstrap consensus tree</li>
                    <li><strong>Tree visualizations</strong> - PNG and PDF plots</li>
                    <li><strong>IQ-TREE log files</strong> - Detailed analysis results</li>
                </ul>
                
                <h3>Recommended Tree Viewers:</h3>
                <ul>
                    <li><a href="https://github.com/rambaut/figtree/">FigTree</a> - Interactive tree visualization</li>
                    <li><a href="https://itol.embl.de/">iTOL</a> - Online tree annotation and display</li>
                    <li><a href="http://dendroscope.org/">Dendroscope</a> - Advanced phylogenetic tree viewer</li>
                    <li>R packages: ape, ggtree, phytools</li>
                </ul>
            </div>
        </div>
        
        <div class="section">
            <h2>Methodology</h2>
            <p>This phylogenetic analysis was performed using IQ-TREE with the following approach:</p>
            <ol>
                <li><strong>Input:</strong> Core genome alignment from pan-genome analysis</li>
                <li><strong>Model Selection:</strong> Automatic model testing to find best-fit evolutionary model</li>
                <li><strong>Tree Construction:</strong> Maximum likelihood inference</li>
                <li><strong>Bootstrap Support:</strong> 1000 ultrafast bootstrap replicates for branch support</li>
                <li><strong>Output:</strong> Publication-ready phylogenetic trees with support values</li>
            </ol>
        </div>
        
        <div class="section">
            <h2>Interpretation</h2>
            <p>The phylogenetic tree shows the evolutionary relationships between your bacterial genomes based on core genome sequences. 
            Bootstrap support values indicate the confidence in each branch, with values ≥70% generally considered well-supported.</p>
            
            {% if stats.n_sites %}
            <p>This analysis used {{ "{:,}".format(stats.n_sites) }} alignment positions from the core genome, providing robust phylogenetic signal for inferring relationships.</p>
            {% endif %}
        </div>
        
        <div class="footer">
            <hr>
            <p><em>Report generated by IQ-TREE Phylogenetic Analysis Pipeline</em></p>
        </div>
    </body>
    </html>
    """
    
    template_data = {
        'date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'stats': stats,
        'tree_files': tree_files
    }
    
    template = Template(html_template)
    html_content = template.render(**template_data)
    
    with open(output_file, 'w') as f:
        f.write(html_content)

def main():
    # Get input and output files from Snakemake
    tree_file = snakemake.input.tree
    consensus_tree = snakemake.input.consensus_tree
    iqtree_log = snakemake.input.iqtree_log
    alignment_info = snakemake.input.alignment_info
    
    tree_plot_output = snakemake.output.tree_plot
    tree_plot_pdf = snakemake.output.tree_plot_pdf
    consensus_plot_output = snakemake.output.consensus_plot
    html_output = snakemake.output.html_report
    tree_stats_output = snakemake.output.tree_stats
    
    print("Creating phylogenetic tree visualizations...")
    
    # Create output directory
    output_dir = Path(html_output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse IQ-TREE statistics
    print("Parsing IQ-TREE results...")
    stats = parse_iqtree_log(iqtree_log)
    
    # Create tree visualizations
    print("Creating maximum likelihood tree plot...")
    success_ml = visualize_tree_matplotlib(tree_file, tree_plot_output, 
                                          "Maximum Likelihood Phylogenetic Tree")
    
    # Create PDF version
    visualize_tree_matplotlib(tree_file, tree_plot_pdf, 
                             "Maximum Likelihood Phylogenetic Tree")
    
    print("Creating consensus tree plot...")
    success_consensus = visualize_tree_matplotlib(consensus_tree, consensus_plot_output,
                                                 "Bootstrap Consensus Tree")
    
    # Create statistics plot
    stats_plot = str(tree_plot_output).replace('.png', '_statistics.png')
    create_tree_statistics_plot(stats, stats_plot)
    
    # Create HTML report
    print("Creating HTML report...")
    tree_files = {
        'tree': Path(tree_file).name,
        'consensus': Path(consensus_tree).name
    }
    
    create_html_report(stats, tree_files, html_output)
    
    # Write tree statistics
    print("Writing tree statistics...")
    with open(tree_stats_output, 'w') as f:
        f.write("Phylogenetic Tree Statistics\n")
        f.write("===========================\n\n")
        
        if 'n_sequences' in stats:
            f.write(f"Number of sequences: {stats['n_sequences']}\n")
        if 'n_sites' in stats:
            f.write(f"Alignment length: {stats['n_sites']:,} bp\n")
        if 'best_model' in stats:
            f.write(f"Best-fit model: {stats['best_model']}\n")
        if 'log_likelihood' in stats:
            f.write(f"Log-likelihood: {stats['log_likelihood']:.2f}\n")
        if 'aic' in stats:
            f.write(f"AIC score: {stats['aic']:.2f}\n")
        if 'bic' in stats:
            f.write(f"BIC score: {stats['bic']:.2f}\n")
        if 'cpu_time' in stats:
            f.write(f"CPU time: {stats['cpu_time']:.1f} seconds\n")
        if 'wall_time' in stats:
            f.write(f"Wall time: {stats['wall_time']:.1f} seconds\n")
        
        f.write(f"\nTree files generated:\n")
        f.write(f"- Maximum likelihood tree: {tree_file}\n")
        f.write(f"- Bootstrap consensus tree: {consensus_tree}\n")
        f.write(f"- Tree visualizations: {tree_plot_output}, {consensus_plot_output}\n")
        f.write(f"- HTML report: {html_output}\n")
    
    print("Phylogenetic analysis visualization completed!")
    
    if success_ml and success_consensus:
        print("Tree visualizations created successfully!")
    else:
        print("Note: Some visualizations may require external tools like FigTree for optimal display")

if __name__ == "__main__":
    main()
