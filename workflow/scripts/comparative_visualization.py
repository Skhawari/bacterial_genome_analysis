#!/usr/bin/env python3
"""
Comparative Genomics Visualization Script
Creates comprehensive reports and visualizations for Roary pan-genome analysis
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
from pathlib import Path
from jinja2 import Template
from collections import Counter

# Set style for matplotlib
try:
    plt.style.use('seaborn-v0_8')
except OSError:
    # Fallback for older matplotlib/seaborn versions
    try:
        plt.style.use('seaborn')
    except OSError:
        plt.style.use('default')
        
sns.set_palette("husl")

def parse_gene_presence_absence(csv_file):
    """Parse Roary gene presence/absence matrix"""
    try:
        df = pd.read_csv(csv_file, sep=',', low_memory=False)
        
        # Get sample columns (exclude annotation columns)
        annotation_cols = ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 
                          'No. sequences', 'Avg sequences per isolate', 'Genome Fragment', 
                          'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']
        
        sample_cols = [col for col in df.columns if col not in annotation_cols]
        
        return df, sample_cols
    except Exception as e:
        print(f"Error parsing gene presence/absence file: {e}")
        return None, []

def analyze_pangenome_structure(df, sample_cols):
    """Analyze pan-genome structure (core, shell, cloud)"""
    n_samples = len(sample_cols)
    
    # Count presence per gene
    gene_counts = []
    for idx, row in df.iterrows():
        presence_count = 0
        for col in sample_cols:
            if pd.notna(row[col]) and row[col] != '':
                presence_count += 1
        gene_counts.append(presence_count)
    
    df['presence_count'] = gene_counts
    
    # Define categories
    core_threshold = 0.99 * n_samples  # 99% of samples
    shell_threshold = 0.15 * n_samples  # 15% of samples
    
    categories = []
    for count in gene_counts:
        if count >= core_threshold:
            categories.append('Core')
        elif count >= shell_threshold:
            categories.append('Shell')
        else:
            categories.append('Cloud')
    
    df['category'] = categories
    
    # Calculate statistics
    stats = {
        'total_genes': len(df),
        'core_genes': sum(1 for cat in categories if cat == 'Core'),
        'shell_genes': sum(1 for cat in categories if cat == 'Shell'),
        'cloud_genes': sum(1 for cat in categories if cat == 'Cloud'),
        'n_samples': n_samples
    }
    
    return df, stats

def create_pangenome_plot(df, stats):
    """Create pan-genome structure visualization"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Pie chart of gene categories
    categories = ['Core', 'Shell', 'Cloud']
    counts = [stats['core_genes'], stats['shell_genes'], stats['cloud_genes']]
    colors = ['#2E86AB', '#A23B72', '#F18F01']
    
    wedges, texts, autotexts = ax1.pie(counts, labels=categories, colors=colors, autopct='%1.1f%%',
                                      startangle=90, textprops={'fontsize': 12})
    ax1.set_title('Pan-genome Structure', fontsize=14, fontweight='bold')
    
    # Add count annotations
    for i, (wedge, count) in enumerate(zip(wedges, counts)):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = 0.7 * np.cos(np.radians(angle))
        y = 0.7 * np.sin(np.radians(angle))
        ax1.annotate(f'{count:,} genes', xy=(x, y), ha='center', va='center',
                    fontsize=10, fontweight='bold', color='white')
    
    # Gene frequency distribution
    frequency_counts = Counter(df['presence_count'])
    frequencies = sorted(frequency_counts.keys())
    counts_freq = [frequency_counts[f] for f in frequencies]
    
    bars = ax2.bar(frequencies, counts_freq, color='#2E86AB', alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Number of Genomes', fontsize=12)
    ax2.set_ylabel('Number of Gene Families', fontsize=12)
    ax2.set_title('Gene Frequency Distribution', fontsize=14, fontweight='bold')
    
    # Add core genome threshold line
    core_threshold = 0.99 * stats['n_samples']
    ax2.axvline(x=core_threshold, color='red', linestyle='--', linewidth=2, 
                label=f'Core threshold (99%)')
    ax2.legend()
    
    # Add value labels on bars for high frequencies
    for bar, freq, count in zip(bars, frequencies, counts_freq):
        if freq >= core_threshold or count > max(counts_freq) * 0.1:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count:,}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig

def create_core_accessory_plot(df, sample_cols, stats):
    """Create core vs accessory genome comparison"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Sample-wise gene counts
    sample_data = []
    for sample in sample_cols:
        total_genes = 0
        core_genes = 0
        shell_genes = 0
        cloud_genes = 0
        
        for idx, row in df.iterrows():
            if pd.notna(row[sample]) and row[sample] != '':
                total_genes += 1
                if row['category'] == 'Core':
                    core_genes += 1
                elif row['category'] == 'Shell':
                    shell_genes += 1
                else:
                    cloud_genes += 1
        
        sample_data.append({
            'Sample': sample,
            'Core': core_genes,
            'Shell': shell_genes,
            'Cloud': cloud_genes,
            'Total': total_genes
        })
    
    df_samples = pd.DataFrame(sample_data)
    
    # Stacked bar chart
    x = np.arange(len(sample_cols))
    width = 0.8
    
    p1 = ax1.bar(x, df_samples['Core'], width, label='Core', color='#2E86AB')
    p2 = ax1.bar(x, df_samples['Shell'], width, bottom=df_samples['Core'], 
                label='Shell', color='#A23B72')
    p3 = ax1.bar(x, df_samples['Cloud'], width, 
                bottom=df_samples['Core'] + df_samples['Shell'],
                label='Cloud', color='#F18F01')
    
    ax1.set_xlabel('Samples', fontsize=12)
    ax1.set_ylabel('Number of Genes', fontsize=12)
    ax1.set_title('Gene Content per Sample', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(sample_cols, rotation=45, ha='right')
    ax1.legend()
    
    # Add total counts on bars
    for i, total in enumerate(df_samples['Total']):
        ax1.text(i, total + 50, f'{total:,}', ha='center', va='bottom', fontsize=9)
    
    # Core genome size variation
    core_genes_per_sample = df_samples['Core']
    ax2.bar(x, core_genes_per_sample, color='#2E86AB', alpha=0.7)
    ax2.set_xlabel('Samples', fontsize=12)
    ax2.set_ylabel('Core Genes Count', fontsize=12)
    ax2.set_title('Core Genome Size Variation', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(sample_cols, rotation=45, ha='right')
    
    # Add mean line
    mean_core = np.mean(core_genes_per_sample)
    ax2.axhline(y=mean_core, color='red', linestyle='--', linewidth=2,
                label=f'Mean: {mean_core:.0f}')
    ax2.legend()
    
    plt.tight_layout()
    return fig

def create_gene_frequency_plot(df, stats):
    """Create detailed gene frequency distribution plot"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create frequency distribution
    frequency_counts = Counter(df['presence_count'])
    frequencies = sorted(frequency_counts.keys())
    counts = [frequency_counts[f] for f in frequencies]
    
    # Create histogram-style plot
    bars = ax.bar(frequencies, counts, alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Color bars by category
    core_threshold = 0.99 * stats['n_samples']
    shell_threshold = 0.15 * stats['n_samples']
    
    for bar, freq in zip(bars, frequencies):
        if freq >= core_threshold:
            bar.set_color('#2E86AB')  # Core = blue
        elif freq >= shell_threshold:
            bar.set_color('#A23B72')  # Shell = purple
        else:
            bar.set_color('#F18F01')  # Cloud = orange
    
    ax.set_xlabel('Number of Genomes Containing Gene Family', fontsize=12)
    ax.set_ylabel('Number of Gene Families', fontsize=12)
    ax.set_title('Gene Family Frequency Distribution', fontsize=14, fontweight='bold')
    
    # Add threshold lines
    ax.axvline(x=core_threshold, color='blue', linestyle='--', linewidth=2,
               label=f'Core threshold (≥{core_threshold:.1f})')
    ax.axvline(x=shell_threshold, color='purple', linestyle='--', linewidth=2,
               label=f'Shell threshold (≥{shell_threshold:.1f})')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2E86AB', label='Core genes'),
        Patch(facecolor='#A23B72', label='Shell genes'),
        Patch(facecolor='#F18F01', label='Cloud genes'),
        ax.lines[0],  # Core threshold line
        ax.lines[1]   # Shell threshold line
    ]
    ax.legend(handles=legend_elements)
    
    plt.tight_layout()
    return fig

def create_html_report(stats, sample_cols, plots_dir):
    """Create comprehensive HTML report for comparative genomics"""
    
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Comparative Genomics Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
            .header { text-align: center; color: #2E86AB; margin-bottom: 30px; }
            .summary-table { width: 100%; border-collapse: collapse; margin: 20px 0; }
            .summary-table th, .summary-table td { border: 1px solid #ddd; padding: 10px; text-align: left; }
            .summary-table th { background-color: #f2f2f2; font-weight: bold; }
            .metric { display: inline-block; margin: 15px; padding: 20px; background: #f8f9fa; border-radius: 8px; }
            .metric-value { font-size: 32px; font-weight: bold; color: #2E86AB; }
            .metric-label { font-size: 14px; color: #666; margin-top: 5px; }
            .plot-container { text-align: center; margin: 30px 0; }
            .section { margin: 40px 0; }
            .core { color: #2E86AB; font-weight: bold; }
            .shell { color: #A23B72; font-weight: bold; }
            .cloud { color: #F18F01; font-weight: bold; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Comparative Genomics Analysis Report</h1>
            <h2>Pan-genome Analysis Results</h2>
            <p>Generated on {{ date }}</p>
        </div>
        
        <div class="section">
            <h2>Pan-genome Overview</h2>
            <div style="text-align: center;">
                <div class="metric">
                    <div class="metric-value">{{ stats.n_samples }}</div>
                    <div class="metric-label">Genomes Analyzed</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{{ "{:,}".format(stats.total_genes) }}</div>
                    <div class="metric-label">Total Gene Families</div>
                </div>
                <div class="metric">
                    <div class="metric-value core">{{ "{:,}".format(stats.core_genes) }}</div>
                    <div class="metric-label">Core Genes (99%+)</div>
                </div>
                <div class="metric">
                    <div class="metric-value shell">{{ "{:,}".format(stats.shell_genes) }}</div>
                    <div class="metric-label">Shell Genes (15-99%)</div>
                </div>
                <div class="metric">
                    <div class="metric-value cloud">{{ "{:,}".format(stats.cloud_genes) }}</div>
                    <div class="metric-label">Cloud Genes (<15%)</div>
                </div>
            </div>
        </div>
        
        <div class="section">
            <h2>Pan-genome Statistics</h2>
            <table class="summary-table">
                <tr>
                    <th>Category</th>
                    <th>Gene Families</th>
                    <th>Percentage</th>
                    <th>Definition</th>
                </tr>
                <tr>
                    <td class="core">Core Genome</td>
                    <td>{{ "{:,}".format(stats.core_genes) }}</td>
                    <td>{{ "%.1f"|format((stats.core_genes / stats.total_genes) * 100) }}%</td>
                    <td>Present in ≥99% of genomes</td>
                </tr>
                <tr>
                    <td class="shell">Shell Genome</td>
                    <td>{{ "{:,}".format(stats.shell_genes) }}</td>
                    <td>{{ "%.1f"|format((stats.shell_genes / stats.total_genes) * 100) }}%</td>
                    <td>Present in 15-99% of genomes</td>
                </tr>
                <tr>
                    <td class="cloud">Cloud Genome</td>
                    <td>{{ "{:,}".format(stats.cloud_genes) }}</td>
                    <td>{{ "%.1f"|format((stats.cloud_genes / stats.total_genes) * 100) }}%</td>
                    <td>Present in <15% of genomes</td>
                </tr>
                <tr style="background-color: #f8f9fa; font-weight: bold;">
                    <td>Total Pan-genome</td>
                    <td>{{ "{:,}".format(stats.total_genes) }}</td>
                    <td>100.0%</td>
                    <td>All gene families identified</td>
                </tr>
            </table>
        </div>
        
        <div class="section">
            <h2>Analysis Summary</h2>
            <p>This pan-genome analysis of <strong>{{ stats.n_samples }} genomes</strong> identified a total of 
            <strong>{{ "{:,}".format(stats.total_genes) }} gene families</strong>. The core genome, consisting of 
            genes present in at least 99% of the genomes, contains <span class="core">{{ "{:,}".format(stats.core_genes) }} gene families</span> 
            ({{ "%.1f"|format((stats.core_genes / stats.total_genes) * 100) }}% of the pan-genome).</p>
            
            <p>The accessory genome is composed of <span class="shell">{{ "{:,}".format(stats.shell_genes) }} shell genes</span> 
            ({{ "%.1f"|format((stats.shell_genes / stats.total_genes) * 100) }}%) and 
            <span class="cloud">{{ "{:,}".format(stats.cloud_genes) }} cloud genes</span> 
            ({{ "%.1f"|format((stats.cloud_genes / stats.total_genes) * 100) }}%), representing the variable 
            genetic content across the analyzed genomes.</p>
        </div>
        
        <div class="section">
            <h2>Ready for Phylogenetic Analysis</h2>
            <p>The core genome alignment has been generated and is ready for phylogenetic tree construction using IQ-TREE. 
            The alignment contains <strong>{{ "{:,}".format(stats.core_genes) }} core gene families</strong> that will provide 
            robust evolutionary signal for inferring relationships between the analyzed genomes.</p>
        </div>
        
        <div class="footer">
            <hr>
            <p><em>Report generated by Comparative Genomics Analysis Pipeline (Roary)</em></p>
        </div>
    </body>
    </html>
    """
    
    template_data = {
        'date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'stats': stats,
        'samples': sample_cols
    }
    
    template = Template(html_template)
    html_content = template.render(**template_data)
    
    return html_content

def main():
    # Get input and output files from Snakemake
    gene_presence_absence = snakemake.input.gene_presence_absence
    html_output = snakemake.output.html_report
    pangenome_plot_output = snakemake.output.pangenome_plot
    core_accessory_plot_output = snakemake.output.core_accessory_plot
    gene_frequency_plot_output = snakemake.output.gene_frequency_plot
    
    print("Parsing gene presence/absence matrix...")
    df, sample_cols = parse_gene_presence_absence(gene_presence_absence)
    
    if df is None:
        print("Error: Could not parse gene presence/absence file")
        return
    
    print(f"Analyzing pan-genome structure for {len(sample_cols)} samples...")
    df, stats = analyze_pangenome_structure(df, sample_cols)
    
    print(f"Pan-genome statistics:")
    print(f"  Total gene families: {stats['total_genes']:,}")
    print(f"  Core genes: {stats['core_genes']:,} ({stats['core_genes']/stats['total_genes']*100:.1f}%)")
    print(f"  Shell genes: {stats['shell_genes']:,} ({stats['shell_genes']/stats['total_genes']*100:.1f}%)")
    print(f"  Cloud genes: {stats['cloud_genes']:,} ({stats['cloud_genes']/stats['total_genes']*100:.1f}%)")
    
    # Create output directory
    output_dir = Path(html_output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    print("Creating pan-genome structure plot...")
    pangenome_fig = create_pangenome_plot(df, stats)
    pangenome_fig.savefig(pangenome_plot_output, dpi=300, bbox_inches='tight')
    plt.close(pangenome_fig)
    
    print("Creating core vs accessory plot...")
    core_accessory_fig = create_core_accessory_plot(df, sample_cols, stats)
    core_accessory_fig.savefig(core_accessory_plot_output, dpi=300, bbox_inches='tight')
    plt.close(core_accessory_fig)
    
    print("Creating gene frequency distribution plot...")
    gene_freq_fig = create_gene_frequency_plot(df, stats)
    gene_freq_fig.savefig(gene_frequency_plot_output, dpi=300, bbox_inches='tight')
    plt.close(gene_freq_fig)
    
    # Generate HTML report
    print("Creating HTML report...")
    html_content = create_html_report(stats, sample_cols, output_dir)
    with open(html_output, 'w') as f:
        f.write(html_content)
    
    print("Comparative genomics visualization completed!")
    print(f"Core genome contains {stats['core_genes']:,} gene families - ready for phylogenetic analysis!")

if __name__ == "__main__":
    main()
