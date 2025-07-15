#!/usr/bin/env python3
"""
Prokka Annotation Visualization Script
Creates comprehensive reports and visualizations for bacterial genome annotations
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
from pathlib import Path
import re
from jinja2 import Template
import weasyprint
from collections import defaultdict, Counter

# Set style for matplotlib
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def parse_prokka_txt(txt_file):
    """Parse Prokka summary statistics from .txt file"""
    stats = {}
    with open(txt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()
                try:
                    stats[key] = int(value)
                except ValueError:
                    stats[key] = value
    return stats

def parse_prokka_gff(gff_file):
    """Extract feature statistics from GFF file"""
    features = defaultdict(int)
    contigs = set()
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('##sequence-region'):
                    contig = line.split()[1]
                    contigs.add(contig)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                feature_type = fields[2]
                features[feature_type] += 1
    
    return dict(features), len(contigs)

def parse_prokka_tsv(tsv_file):
    """Parse functional annotation data from TSV file"""
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        
        # Count functional categories
        functional_stats = {
            'with_gene_name': df['gene'].notna().sum(),
            'with_ec_number': df['EC_number'].notna().sum(),
            'with_cog': df['COG'].notna().sum(),
            'hypothetical': df['product'].str.contains('hypothetical', case=False, na=False).sum(),
            'total_features': len(df)
        }
        
        # Extract COG categories
        cog_categories = []
        for cog in df['COG'].dropna():
            if isinstance(cog, str) and cog.startswith('COG'):
                cog_categories.append(cog)
        
        return functional_stats, cog_categories
    except Exception as e:
        print(f"Error parsing {tsv_file}: {e}")
        return {}, []

def create_summary_plots(all_stats, sample_names):
    """Create summary visualization plots"""
    
    # Prepare data for plotting
    metrics = ['contigs', 'bases', 'CDS', 'gene', 'rRNA', 'tRNA']
    plot_data = []
    
    for sample in sample_names:
        stats = all_stats[sample]['txt_stats']
        for metric in metrics:
            if metric in stats:
                plot_data.append({
                    'Sample': sample,
                    'Metric': metric,
                    'Value': stats[metric]
                })
    
    df_plot = pd.DataFrame(plot_data)
    
    # Create subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Prokka Annotation Summary - All Samples', fontsize=16, fontweight='bold')
    
    for i, metric in enumerate(metrics):
        row = i // 3
        col = i % 3
        ax = axes[row, col]
        
        metric_data = df_plot[df_plot['Metric'] == metric]
        if not metric_data.empty:
            bars = ax.bar(metric_data['Sample'], metric_data['Value'], 
                         color=sns.color_palette("husl", len(sample_names)))
            ax.set_title(f'{metric.upper()} Count', fontweight='bold')
            ax.set_xlabel('Sample')
            ax.set_ylabel('Count')
            ax.tick_params(axis='x', rotation=45)
            
            # Add value labels on bars
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height):,}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig

def create_functional_plot(all_stats, sample_names):
    """Create functional annotation comparison plot"""
    
    # Prepare functional annotation data
    func_data = []
    for sample in sample_names:
        func_stats = all_stats[sample]['functional_stats']
        total = func_stats.get('total_features', 1)
        
        func_data.append({
            'Sample': sample,
            'With Gene Name': func_stats.get('with_gene_name', 0),
            'With EC Number': func_stats.get('with_ec_number', 0), 
            'With COG': func_stats.get('with_cog', 0),
            'Hypothetical': func_stats.get('hypothetical', 0),
            'Total Features': total
        })
    
    df_func = pd.DataFrame(func_data)
    
    # Calculate percentages
    for col in ['With Gene Name', 'With EC Number', 'With COG', 'Hypothetical']:
        df_func[f'{col} (%)'] = (df_func[col] / df_func['Total Features']) * 100
    
    # Create grouped bar plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Absolute counts
    categories = ['With Gene Name', 'With EC Number', 'With COG', 'Hypothetical']
    x = np.arange(len(sample_names))
    width = 0.2
    
    for i, category in enumerate(categories):
        offset = (i - 1.5) * width
        bars = ax1.bar(x + offset, df_func[category], width, label=category)
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    ax1.set_title('Functional Annotation Counts', fontweight='bold')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Number of Features')
    ax1.set_xticks(x)
    ax1.set_xticklabels(sample_names, rotation=45)
    ax1.legend()
    
    # Percentages
    for i, category in enumerate(categories):
        offset = (i - 1.5) * width
        bars = ax2.bar(x + offset, df_func[f'{category} (%)'], width, label=category)
        
        # Add percentage labels
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.1f}%', ha='center', va='bottom', fontsize=8)
    
    ax2.set_title('Functional Annotation Percentages', fontweight='bold')
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Percentage of Features')
    ax2.set_xticks(x)
    ax2.set_xticklabels(sample_names, rotation=45)
    ax2.legend()
    
    plt.tight_layout()
    return fig

def create_interactive_plotly_report(all_stats, sample_names):
    """Create interactive Plotly visualizations"""
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=['Genome Size Distribution', 'Gene Count Comparison', 
                       'Feature Type Distribution', 'Functional Annotation Quality'],
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"type": "pie"}, {"secondary_y": False}]]
    )
    
    # 1. Genome size scatter plot
    genome_sizes = [all_stats[sample]['txt_stats'].get('bases', 0) for sample in sample_names]
    gene_counts = [all_stats[sample]['txt_stats'].get('CDS', 0) for sample in sample_names]
    
    fig.add_trace(
        go.Scatter(x=sample_names, y=genome_sizes, mode='markers+lines',
                  name='Genome Size (bp)', marker=dict(size=12)),
        row=1, col=1
    )
    
    # 2. Gene count comparison
    metrics = ['CDS', 'rRNA', 'tRNA']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, metric in enumerate(metrics):
        values = [all_stats[sample]['txt_stats'].get(metric, 0) for sample in sample_names]
        fig.add_trace(
            go.Bar(x=sample_names, y=values, name=metric, marker_color=colors[i]),
            row=1, col=2
        )
    
    # 3. Feature type pie chart (using first sample as example)
    sample1_features = all_stats[sample_names[0]]['gff_features']
    fig.add_trace(
        go.Pie(labels=list(sample1_features.keys()), 
               values=list(sample1_features.values()),
               name=f"Features - {sample_names[0]}"),
        row=2, col=1
    )
    
    # 4. Functional annotation quality
    func_categories = ['With Gene Name', 'With EC Number', 'With COG', 'Hypothetical']
    
    for category in func_categories:
        values = []
        for sample in sample_names:
            func_stats = all_stats[sample]['functional_stats']
            total = func_stats.get('total_features', 1)
            if 'Gene Name' in category:
                val = func_stats.get('with_gene_name', 0)
            elif 'EC Number' in category:
                val = func_stats.get('with_ec_number', 0)
            elif 'COG' in category:
                val = func_stats.get('with_cog', 0)
            else:
                val = func_stats.get('hypothetical', 0)
            
            values.append((val / total) * 100 if total > 0 else 0)
        
        fig.add_trace(
            go.Bar(x=sample_names, y=values, name=category),
            row=2, col=2
        )
    
    # Update layout
    fig.update_layout(
        height=800,
        title_text="Prokka Annotation Analysis - Interactive Dashboard",
        title_x=0.5,
        showlegend=True
    )
    
    # Update axes labels
    fig.update_xaxes(title_text="Sample", row=1, col=1)
    fig.update_yaxes(title_text="Genome Size (bp)", row=1, col=1)
    fig.update_xaxes(title_text="Sample", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=2)
    fig.update_xaxes(title_text="Sample", row=2, col=2)
    fig.update_yaxes(title_text="Percentage (%)", row=2, col=2)
    
    return fig

def create_html_report(all_stats, sample_names, plots_dir, summary_plot_path, functional_plot_path):
    """Create comprehensive HTML report with embedded plots"""
    
    # Convert plot images to base64 for embedding
    import base64
    
    def image_to_base64(image_path):
        try:
            with open(image_path, "rb") as img_file:
                return base64.b64encode(img_file.read()).decode()
        except:
            return None
    
    summary_plot_b64 = image_to_base64(summary_plot_path)
    functional_plot_b64 = image_to_base64(functional_plot_path)
    
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Prokka Annotation Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
            .header { text-align: center; color: #2E86AB; margin-bottom: 30px; }
            .summary-table { width: 100%; border-collapse: collapse; margin: 20px 0; }
            .summary-table th, .summary-table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            .summary-table th { background-color: #f2f2f2; font-weight: bold; }
            .sample-section { margin: 30px 0; padding: 20px; border: 1px solid #eee; border-radius: 8px; }
            .metric { display: inline-block; margin: 10px; padding: 15px; background: #f8f9fa; border-radius: 5px; }
            .metric-value { font-size: 24px; font-weight: bold; color: #2E86AB; }
            .metric-label { font-size: 14px; color: #666; }
            .plot-container { text-align: center; margin: 30px 0; page-break-inside: avoid; }
            .plot-image { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 8px; }
            .plot-title { font-size: 18px; font-weight: bold; margin-bottom: 15px; color: #333; }
            @media print {
                .plot-container { page-break-inside: avoid; }
                body { margin: 20px; }
            }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Prokka Genome Annotation Analysis Report</h1>
            <h2>Klebsiella pneumoniae Strains</h2>
            <p>Generated on {{ date }}</p>
        </div>
        
        <h2>Executive Summary</h2>
        <table class="summary-table">
            <tr>
                <th>Sample</th>
                <th>Contigs</th>
                <th>Genome Size (Mb)</th>
                <th>Total Genes</th>
                <th>Protein-coding</th>
                <th>rRNA</th>
                <th>tRNA</th>
                <th>Functional (%)</th>
            </tr>
            {% for sample in samples %}
            <tr>
                <td><strong>{{ sample }}</strong></td>
                <td>{{ stats[sample].contigs }}</td>
                <td>{{ "%.2f"|format(stats[sample].bases / 1000000) }}</td>
                <td>{{ stats[sample].gene }}</td>
                <td>{{ stats[sample].CDS }}</td>
                <td>{{ stats[sample].rRNA }}</td>
                <td>{{ stats[sample].tRNA }}</td>
                <td>{{ "%.1f"|format(stats[sample].functional_percent) }}%</td>
            </tr>
            {% endfor %}
        </table>
        
        <h2>Detailed Analysis</h2>
        
        {% for sample in samples %}
        <div class="sample-section">
            <h3>Sample {{ sample }}</h3>
            <div class="metrics-container">
                <div class="metric">
                    <div class="metric-value">{{ stats[sample].contigs }}</div>
                    <div class="metric-label">Contigs</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{{ "%.2f"|format(stats[sample].bases / 1000000) }}</div>
                    <div class="metric-label">Genome Size (Mb)</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{{ stats[sample].CDS }}</div>
                    <div class="metric-label">Protein-coding Genes</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{{ stats[sample].rRNA + stats[sample].tRNA }}</div>
                    <div class="metric-label">RNA Genes</div>
                </div>
            </div>
        </div>
        {% endfor %}
        
        <div class="plot-container">
            <div class="plot-title">Annotation Summary - All Samples</div>
            {% if summary_plot_b64 %}
            <img src="data:image/png;base64,{{ summary_plot_b64 }}" class="plot-image" alt="Annotation Summary Plot">
            {% else %}
            <p>Summary plot not available</p>
            {% endif %}
        </div>
        
        <div class="plot-container">
            <div class="plot-title">Functional Annotation Analysis</div>
            {% if functional_plot_b64 %}
            <img src="data:image/png;base64,{{ functional_plot_b64 }}" class="plot-image" alt="Functional Categories Plot">
            {% else %}
            <p>Functional plot not available</p>
            {% endif %}
        </div>
        
        <div class="footer">
            <hr>
            <p><em>Report generated by Prokka Annotation Analysis Pipeline</em></p>
        </div>
    </body>
    </html>
    """
    
    # Prepare data for template
    template_data = {
        'date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'samples': sample_names,
        'stats': {},
        'summary_plot_b64': summary_plot_b64,
        'functional_plot_b64': functional_plot_b64
    }
    
    for sample in sample_names:
        sample_stats = all_stats[sample]['txt_stats']
        func_stats = all_stats[sample]['functional_stats']
        
        total_features = func_stats.get('total_features', 1)
        functional_count = (func_stats.get('with_gene_name', 0) + 
                          func_stats.get('with_ec_number', 0) + 
                          func_stats.get('with_cog', 0))
        functional_percent = (functional_count / total_features) * 100 if total_features > 0 else 0
        
        template_data['stats'][sample] = {
            'contigs': sample_stats.get('contigs', 0),
            'bases': sample_stats.get('bases', 0),
            'gene': sample_stats.get('gene', 0),
            'CDS': sample_stats.get('CDS', 0),
            'rRNA': sample_stats.get('rRNA', 0),
            'tRNA': sample_stats.get('tRNA', 0),
            'functional_percent': functional_percent
        }
    
    # Render template
    template = Template(html_template)
    html_content = template.render(**template_data)
    
    return html_content

def main():
    # Get input and output files from Snakemake
    gff_files = snakemake.input.gff_files
    txt_files = snakemake.input.txt_files
    tsv_files = snakemake.input.tsv_files
    
    html_output = snakemake.output.html_report
    pdf_output = snakemake.output.pdf_report
    summary_plot_output = snakemake.output.summary_plot
    functional_plot_output = snakemake.output.functional_plot
    
    # Extract sample names
    sample_names = []
    for txt_file in txt_files:
        sample_name = Path(txt_file).stem
        sample_names.append(sample_name)
    
    print(f"Processing {len(sample_names)} samples: {sample_names}")
    
    # Parse all data
    all_stats = {}
    
    for i, sample in enumerate(sample_names):
        print(f"Processing sample {sample}...")
        
        # Parse files
        txt_stats = parse_prokka_txt(txt_files[i])
        gff_features, contig_count = parse_prokka_gff(gff_files[i])
        functional_stats, cog_categories = parse_prokka_tsv(tsv_files[i])
        
        all_stats[sample] = {
            'txt_stats': txt_stats,
            'gff_features': gff_features,
            'functional_stats': functional_stats,
            'cog_categories': cog_categories
        }
    
    # Create output directory
    output_dir = Path(html_output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plots first
    print("Creating summary plot...")
    summary_fig = create_summary_plots(all_stats, sample_names)
    summary_fig.savefig(summary_plot_output, dpi=300, bbox_inches='tight')
    plt.close(summary_fig)
    
    print("Creating functional annotation plot...")
    functional_fig = create_functional_plot(all_stats, sample_names)
    functional_fig.savefig(functional_plot_output, dpi=300, bbox_inches='tight')
    plt.close(functional_fig)
    
    # Generate HTML report with embedded plots
    print("Creating HTML report with embedded plots...")
    html_content = create_html_report(all_stats, sample_names, output_dir, 
                                    summary_plot_output, functional_plot_output)
    with open(html_output, 'w') as f:
        f.write(html_content)
    
    # Generate PDF report
    print("Creating PDF report...")
    try:
        weasyprint.HTML(string=html_content).write_pdf(pdf_output)
        print(f"PDF report saved to: {pdf_output}")
    except Exception as e:
        print(f"Warning: Could not generate PDF: {e}")
        # Create a simple PDF fallback
        with open(pdf_output.replace('.pdf', '_fallback.txt'), 'w') as f:
            f.write("PDF generation failed. Please check HTML report.\n")
            f.write(f"Error: {e}\n")
    
    # Generate interactive Plotly report
    print("Creating interactive visualizations...")
    try:
        plotly_fig = create_interactive_plotly_report(all_stats, sample_names)
        plotly_html = html_output.replace('.html', '_interactive.html')
        plotly_fig.write_html(plotly_html)
        print(f"Interactive report saved to: {plotly_html}")
    except Exception as e:
        print(f"Warning: Could not generate interactive plots: {e}")
    
    print("Visualization report generation completed!")

if __name__ == "__main__":
    main()
