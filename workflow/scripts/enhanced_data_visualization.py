#!/usr/bin/env python3
"""
Enhanced Data Visualization Script
Creates comprehensive statistical plots including boxplots, violin plots, and distribution plots
for bacterial genome analysis data
"""

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import json
    from scipy import stats
    import warnings
except ImportError as e:
    print(f"ERROR: Missing required Python package: {e}")
    print("Please ensure all visualization dependencies are installed.")
    print("This script requires: matplotlib, seaborn, pandas, numpy, scipy")
    exit(1)

warnings.filterwarnings('ignore')

# Set style for publication-quality plots
try:
    if 'seaborn-v0_8' in plt.style.available:
        plt.style.use('seaborn-v0_8')
    elif 'seaborn' in plt.style.available:
        plt.style.use('seaborn')
    else:
        plt.style.use('default')
except:
    plt.style.use('default')

sns.set_palette("husl")

def load_assembly_summary_stats(assembly_summary_file):
    """Load assembly summary statistics for boxplot analysis"""
    assembly_data = []
    
    try:
        with open(assembly_summary_file, 'r') as f:
            lines = f.readlines()
            
        # Parse the assembly summary file
        current_sample = None
        for line in lines:
            line = line.strip()
            if line.startswith('Sample:'):
                current_sample = line.split(':')[1].strip()
                sample_data = {'Sample': current_sample}
            elif current_sample and ':' in line:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()
                
                # Extract numeric values
                try:
                    if 'N50' in key:
                        sample_data['N50'] = int(value.replace(',', ''))
                    elif 'Total length' in key:
                        sample_data['Total_length'] = int(value.replace(',', ''))
                    elif 'Number of contigs' in key:
                        sample_data['Contigs'] = int(value.replace(',', ''))
                    elif 'Largest contig' in key:
                        sample_data['Largest_contig'] = int(value.replace(',', ''))
                    elif 'GC content' in key:
                        sample_data['GC_content'] = float(value.replace('%', ''))
                except ValueError:
                    pass
            elif line == '' and current_sample:
                # End of sample data
                assembly_data.append(sample_data)
                current_sample = None
        
        # Add last sample if file doesn't end with empty line
        if current_sample:
            assembly_data.append(sample_data)
            
    except FileNotFoundError:
        print(f"Warning: Assembly summary file not found: {assembly_summary_file}")
        
    return pd.DataFrame(assembly_data)

def load_busco_stats(samples):
    """Load BUSCO completeness statistics"""
    busco_data = []
    
    for sample in samples:
        busco_file = f"results/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.txt"
        try:
            with open(busco_file, 'r') as f:
                content = f.read()
                stats_dict = {'Sample': sample}
                
                # Parse BUSCO percentages
                for line in content.split('\n'):
                    if 'Complete BUSCOs' in line and '%' in line:
                        complete = float(line.split('(')[1].split('%')[0])
                        stats_dict['Complete_BUSCOs'] = complete
                    elif 'Complete and single-copy BUSCOs' in line and '%' in line:
                        single_copy = float(line.split('(')[1].split('%')[0])
                        stats_dict['Single_copy_BUSCOs'] = single_copy
                    elif 'Complete and duplicated BUSCOs' in line and '%' in line:
                        duplicated = float(line.split('(')[1].split('%')[0])
                        stats_dict['Duplicated_BUSCOs'] = duplicated
                    elif 'Fragmented BUSCOs' in line and '%' in line:
                        fragmented = float(line.split('(')[1].split('%')[0])
                        stats_dict['Fragmented_BUSCOs'] = fragmented
                    elif 'Missing BUSCOs' in line and '%' in line:
                        missing = float(line.split('(')[1].split('%')[0])
                        stats_dict['Missing_BUSCOs'] = missing
                
                busco_data.append(stats_dict)
                
        except FileNotFoundError:
            print(f"Warning: BUSCO report not found for {sample}")
            continue
    
    return pd.DataFrame(busco_data)

def load_annotation_stats(samples):
    """Load annotation statistics from Prokka GFF files"""
    annotation_data = []
    
    for sample in samples:
        try:
            stats_dict = {'Sample': sample}
            
            # Get stats from GFF file
            gff_file = f"results/annotation/{sample}/{sample}.gff"
            if Path(gff_file).exists():
                with open(gff_file, 'r') as f:
                    lines = f.readlines()
                    
                    # Count different feature types
                    cds_count = len([l for l in lines if '\tCDS\t' in l])
                    gene_count = len([l for l in lines if '\tgene\t' in l])
                    trna_count = len([l for l in lines if '\ttRNA\t' in l])
                    rrna_count = len([l for l in lines if '\trRNA\t' in l])
                    
                    stats_dict['CDS'] = cds_count
                    stats_dict['Genes'] = gene_count
                    stats_dict['tRNA'] = trna_count
                    stats_dict['rRNA'] = rrna_count
            
            # Try to get additional stats from Prokka log file if available
            prokka_log = f"logs/prokka_{sample}.log"
            if Path(prokka_log).exists():
                try:
                    with open(prokka_log, 'r') as f:
                        content = f.read()
                        
                        # Extract annotation counts from log
                        for line in content.split('\n'):
                            if 'CDS:' in line:
                                try:
                                    cds_count = int(line.split('CDS:')[1].strip())
                                    stats_dict['CDS'] = cds_count
                                except:
                                    pass
                            elif 'tRNA:' in line:
                                try:
                                    trna_count = int(line.split('tRNA:')[1].strip())
                                    stats_dict['tRNA'] = trna_count
                                except:
                                    pass
                            elif 'rRNA:' in line:
                                try:
                                    rrna_count = int(line.split('rRNA:')[1].strip())
                                    stats_dict['rRNA'] = rrna_count
                                except:
                                    pass
                except FileNotFoundError:
                    pass  # Log file not available, use GFF counts
            
            annotation_data.append(stats_dict)
            
        except Exception as e:
            print(f"Warning: Could not load annotation stats for {sample}: {e}")
            continue
    
    return pd.DataFrame(annotation_data)

def create_assembly_quality_plots(df_assembly, df_busco):
    """Create comprehensive assembly quality plots"""
    fig = plt.figure(figsize=(20, 16))
    
    # Create a 3x3 grid for multiple plots
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. N50 Distribution - Boxplot and violin plot
    ax1 = fig.add_subplot(gs[0, 0])
    if 'N50' in df_assembly.columns and not df_assembly['N50'].isna().all():
        sns.boxplot(y=df_assembly['N50'], ax=ax1, color='lightblue')
        ax1.set_title('N50 Distribution', fontsize=12, fontweight='bold')
        ax1.set_ylabel('N50 (bp)')
        ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x/1000)}K'))
    else:
        ax1.text(0.5, 0.5, 'N50 data\nnot available', ha='center', va='center', transform=ax1.transAxes)
        ax1.set_title('N50 Distribution', fontsize=12, fontweight='bold')
    
    # 2. Genome Size Distribution - Violin plot
    ax2 = fig.add_subplot(gs[0, 1])
    if 'Total_length' in df_assembly.columns and not df_assembly['Total_length'].isna().all():
        sns.violinplot(y=df_assembly['Total_length'], ax=ax2, color='lightgreen')
        ax2.set_title('Genome Size Distribution', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Total Length (bp)')
        ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f}M'))
    else:
        ax2.text(0.5, 0.5, 'Genome size data\nnot available', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Genome Size Distribution', fontsize=12, fontweight='bold')
    
    # 3. Contig Count Distribution - Boxplot
    ax3 = fig.add_subplot(gs[0, 2])
    if 'Contigs' in df_assembly.columns and not df_assembly['Contigs'].isna().all():
        sns.boxplot(y=df_assembly['Contigs'], ax=ax3, color='lightcoral')
        ax3.set_title('Contig Count Distribution', fontsize=12, fontweight='bold')
        ax3.set_ylabel('Number of Contigs')
    else:
        ax3.text(0.5, 0.5, 'Contig count data\nnot available', ha='center', va='center', transform=ax3.transAxes)
        ax3.set_title('Contig Count Distribution', fontsize=12, fontweight='bold')
    
    # 4. GC Content Distribution - Histogram with KDE
    ax4 = fig.add_subplot(gs[1, 0])
    if 'GC_content' in df_assembly.columns and not df_assembly['GC_content'].isna().all():
        sns.histplot(df_assembly['GC_content'], kde=True, ax=ax4, color='gold', alpha=0.7)
        ax4.set_title('GC Content Distribution', fontsize=12, fontweight='bold')
        ax4.set_xlabel('GC Content (%)')
        ax4.set_ylabel('Frequency')
    else:
        ax4.text(0.5, 0.5, 'GC content data\nnot available', ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('GC Content Distribution', fontsize=12, fontweight='bold')
    
    # 5. BUSCO Completeness - Multiple boxplots
    ax5 = fig.add_subplot(gs[1, 1])
    if not df_busco.empty and 'Complete_BUSCOs' in df_busco.columns:
        busco_metrics = ['Complete_BUSCOs', 'Single_copy_BUSCOs', 'Duplicated_BUSCOs', 
                        'Fragmented_BUSCOs', 'Missing_BUSCOs']
        busco_data = []
        for metric in busco_metrics:
            if metric in df_busco.columns:
                for value in df_busco[metric].dropna():
                    busco_data.append({'Metric': metric.replace('_BUSCOs', ''), 'Percentage': value})
        
        if busco_data:
            df_busco_plot = pd.DataFrame(busco_data)
            sns.boxplot(data=df_busco_plot, x='Metric', y='Percentage', ax=ax5)
            ax5.set_title('BUSCO Completeness Distribution', fontsize=12, fontweight='bold')
            ax5.set_ylabel('Percentage (%)')
            ax5.tick_params(axis='x', rotation=45)
    else:
        ax5.text(0.5, 0.5, 'BUSCO data\nnot available', ha='center', va='center', transform=ax5.transAxes)
        ax5.set_title('BUSCO Completeness Distribution', fontsize=12, fontweight='bold')
    
    # 6. Assembly Quality vs Completeness - Scatter plot
    ax6 = fig.add_subplot(gs[1, 2])
    if ('N50' in df_assembly.columns and 'Complete_BUSCOs' in df_busco.columns and 
        not df_assembly['N50'].isna().all() and not df_busco['Complete_BUSCOs'].isna().all()):
        
        # Merge dataframes on sample name
        merged_df = pd.merge(df_assembly[['Sample', 'N50']], 
                           df_busco[['Sample', 'Complete_BUSCOs']], on='Sample', how='inner')
        
        if not merged_df.empty:
            scatter = ax6.scatter(merged_df['N50'], merged_df['Complete_BUSCOs'], 
                                alpha=0.7, s=100, c=range(len(merged_df)), cmap='viridis')
            ax6.set_xlabel('N50 (bp)')
            ax6.set_ylabel('BUSCO Completeness (%)')
            ax6.set_title('Assembly Quality vs Completeness', fontsize=12, fontweight='bold')
            ax6.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x/1000)}K'))
            
            # Add sample labels
            for i, row in merged_df.iterrows():
                ax6.annotate(row['Sample'], (row['N50'], row['Complete_BUSCOs']), 
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
    else:
        ax6.text(0.5, 0.5, 'Insufficient data for\nquality comparison', 
                ha='center', va='center', transform=ax6.transAxes)
        ax6.set_title('Assembly Quality vs Completeness', fontsize=12, fontweight='bold')
    
    # 7. Sample Comparison - Multiple metrics
    ax7 = fig.add_subplot(gs[2, :])
    if not df_assembly.empty and 'Sample' in df_assembly.columns:
        # Prepare data for comparison
        comparison_data = []
        metrics_to_plot = ['N50', 'Total_length', 'Contigs', 'GC_content']
        
        for _, row in df_assembly.iterrows():
            sample = row['Sample']
            for metric in metrics_to_plot:
                if metric in row and pd.notna(row[metric]):
                    # Normalize values for better comparison
                    if metric == 'N50':
                        value = row[metric] / 1000  # Convert to KB
                        unit = 'KB'
                    elif metric == 'Total_length':
                        value = row[metric] / 1e6  # Convert to MB
                        unit = 'MB'
                    elif metric == 'Contigs':
                        value = row[metric]
                        unit = 'Count'
                    elif metric == 'GC_content':
                        value = row[metric]
                        unit = '%'
                    
                    comparison_data.append({
                        'Sample': sample,
                        'Metric': f'{metric.replace("_", " ")} ({unit})',
                        'Value': value
                    })
        
        if comparison_data:
            df_comparison = pd.DataFrame(comparison_data)
            sns.boxplot(data=df_comparison, x='Metric', y='Value', ax=ax7)
            ax7.set_title('Assembly Metrics Distribution Across Samples', fontsize=14, fontweight='bold')
            ax7.tick_params(axis='x', rotation=45)
            ax7.set_ylabel('Value')
    else:
        ax7.text(0.5, 0.5, 'Assembly data not available', ha='center', va='center', transform=ax7.transAxes)
        ax7.set_title('Assembly Metrics Distribution Across Samples', fontsize=14, fontweight='bold')
    
    plt.suptitle('Comprehensive Assembly Quality Analysis', fontsize=16, fontweight='bold')
    return fig

def create_annotation_statistics_plots(df_annotation):
    """Create annotation statistics plots"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Genome Annotation Statistical Analysis', fontsize=16, fontweight='bold')
    
    # 1. Gene Count Distribution - Violin plot
    ax1 = axes[0, 0]
    if 'CDS' in df_annotation.columns and not df_annotation['CDS'].isna().all():
        sns.violinplot(y=df_annotation['CDS'], ax=ax1, color='lightblue')
        ax1.set_title('CDS Count Distribution', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Number of CDS')
    else:
        ax1.text(0.5, 0.5, 'CDS data\nnot available', ha='center', va='center', transform=ax1.transAxes)
        ax1.set_title('CDS Count Distribution', fontsize=12, fontweight='bold')
    
    # 2. RNA Features - Boxplot comparison
    ax2 = axes[0, 1]
    rna_data = []
    rna_types = ['tRNA', 'rRNA']
    for rna_type in rna_types:
        if rna_type in df_annotation.columns:
            for value in df_annotation[rna_type].dropna():
                rna_data.append({'RNA_Type': rna_type, 'Count': value})
    
    if rna_data:
        df_rna = pd.DataFrame(rna_data)
        sns.boxplot(data=df_rna, x='RNA_Type', y='Count', ax=ax2)
        ax2.set_title('RNA Features Distribution', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Count')
    else:
        ax2.text(0.5, 0.5, 'RNA data\nnot available', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('RNA Features Distribution', fontsize=12, fontweight='bold')
    
    # 3. Sample Comparison - All annotation features
    ax3 = axes[1, 0]
    if not df_annotation.empty:
        annotation_metrics = ['CDS', 'tRNA', 'rRNA', 'Genes']
        sample_comparison = []
        
        for _, row in df_annotation.iterrows():
            sample = row['Sample']
            for metric in annotation_metrics:
                if metric in row and pd.notna(row[metric]):
                    sample_comparison.append({
                        'Sample': sample,
                        'Feature_Type': metric,
                        'Count': row[metric]
                    })
        
        if sample_comparison:
            df_sample_comp = pd.DataFrame(sample_comparison)
            sns.barplot(data=df_sample_comp, x='Sample', y='Count', hue='Feature_Type', ax=ax3)
            ax3.set_title('Annotation Features by Sample', fontsize=12, fontweight='bold')
            ax3.tick_params(axis='x', rotation=45)
            ax3.legend(title='Feature Type')
    else:
        ax3.text(0.5, 0.5, 'Annotation data\nnot available', ha='center', va='center', transform=ax3.transAxes)
        ax3.set_title('Annotation Features by Sample', fontsize=12, fontweight='bold')
    
    # 4. Gene Density Analysis
    ax4 = axes[1, 1]
    # This would require assembly length data - placeholder for now
    ax4.text(0.5, 0.5, 'Gene density analysis\n(requires assembly data)', 
            ha='center', va='center', transform=ax4.transAxes)
    ax4.set_title('Gene Density Analysis', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    return fig

def create_summary_statistics_report(df_assembly, df_busco, df_annotation, output_file):
    """Create comprehensive summary statistics report"""
    
    # Calculate descriptive statistics
    stats_summary = {}
    
    # Assembly statistics
    if not df_assembly.empty:
        assembly_stats = {}
        numeric_cols = df_assembly.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            if col != 'Sample':
                data = df_assembly[col].dropna()
                if len(data) > 0:
                    assembly_stats[col] = {
                        'mean': float(data.mean()),
                        'median': float(data.median()),
                        'std': float(data.std()),
                        'min': float(data.min()),
                        'max': float(data.max()),
                        'count': int(len(data))
                    }
        stats_summary['assembly'] = assembly_stats
    
    # BUSCO statistics
    if not df_busco.empty:
        busco_stats = {}
        numeric_cols = df_busco.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            if col != 'Sample':
                data = df_busco[col].dropna()
                if len(data) > 0:
                    busco_stats[col] = {
                        'mean': float(data.mean()),
                        'median': float(data.median()),
                        'std': float(data.std()),
                        'min': float(data.min()),
                        'max': float(data.max()),
                        'count': int(len(data))
                    }
        stats_summary['busco'] = busco_stats
    
    # Annotation statistics
    if not df_annotation.empty:
        annotation_stats = {}
        numeric_cols = df_annotation.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            if col != 'Sample':
                data = df_annotation[col].dropna()
                if len(data) > 0:
                    annotation_stats[col] = {
                        'mean': float(data.mean()),
                        'median': float(data.median()),
                        'std': float(data.std()),
                        'min': float(data.min()),
                        'max': float(data.max()),
                        'count': int(len(data))
                    }
        stats_summary['annotation'] = annotation_stats
    
    # Save to JSON
    with open(output_file, 'w') as f:
        json.dump(stats_summary, f, indent=2)
    
    return stats_summary

def create_comprehensive_pdf_report(df_assembly, df_busco, df_annotation, assembly_summary_df, output_file):
    """Create a comprehensive PDF report with all statistical analysis"""
    from matplotlib.backends.backend_pdf import PdfPages
    
    with PdfPages(output_file) as pdf:
        # Page 1: Assembly Quality Analysis
        fig_assembly = create_assembly_quality_plots(df_assembly, df_busco)
        pdf.savefig(fig_assembly, bbox_inches='tight')
        plt.close(fig_assembly)
        
        # Page 2: Assembly Summary Boxplots
        if not assembly_summary_df.empty:
            fig_summary = create_assembly_summary_boxplots(assembly_summary_df)
            pdf.savefig(fig_summary, bbox_inches='tight')
            plt.close(fig_summary)
        
        # Page 3: Annotation Statistics
        fig_annotation = create_annotation_statistics_plots(df_annotation)
        pdf.savefig(fig_annotation, bbox_inches='tight')
        plt.close(fig_annotation)
        
        # Page 4: Summary Tables and Statistics
        fig_tables = create_summary_tables_figure(df_assembly, df_busco, df_annotation, assembly_summary_df)
        pdf.savefig(fig_tables, bbox_inches='tight')
        plt.close(fig_tables)

def create_assembly_summary_boxplots(df_summary):
    """Create boxplots for assembly summary statistics"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Assembly Summary Statistical Analysis', fontsize=16, fontweight='bold')
    
    # Define metrics to plot
    metrics = [
        ('N50', 'N50 (bp)', lambda x: x/1000, 'KB'),
        ('Total_length', 'Genome Size (bp)', lambda x: x/1e6, 'MB'),
        ('Contigs', 'Number of Contigs', lambda x: x, 'Count'),
        ('Largest_contig', 'Largest Contig (bp)', lambda x: x/1000, 'KB'),
        ('GC_content', 'GC Content (%)', lambda x: x, '%')
    ]
    
    for i, (metric, title, transform, unit) in enumerate(metrics):
        if i >= 6:  # Only 6 subplots available
            break
            
        row = i // 3
        col = i % 3
        ax = axes[row, col]
        
        if metric in df_summary.columns and not df_summary[metric].isna().all():
            data = df_summary[metric].dropna()
            transformed_data = data.apply(transform)
            
            # Create boxplot
            box_plot = ax.boxplot(transformed_data, patch_artist=True)
            box_plot['boxes'][0].set_facecolor('lightblue')
            
            # Add individual points
            y_pos = np.random.normal(1, 0.04, len(transformed_data))
            ax.scatter(y_pos, transformed_data, alpha=0.6, color='red', s=30)
            
            # Add statistics text
            stats_text = f"""
            Mean: {transformed_data.mean():.1f} {unit}
            Median: {transformed_data.median():.1f} {unit}
            Std: {transformed_data.std():.1f} {unit}
            Min: {transformed_data.min():.1f} {unit}
            Max: {transformed_data.max():.1f} {unit}
            """
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', fontsize=8,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            ax.set_title(title, fontweight='bold')
            ax.set_ylabel(f'{unit}')
            ax.set_xticklabels(['All Samples'])
        else:
            ax.text(0.5, 0.5, f'{title}\nData not available', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title, fontweight='bold')
    
    # Remove empty subplot
    if len(metrics) < 6:
        axes[1, 2].remove()
    
    plt.tight_layout()
    return fig

def create_summary_tables_figure(df_assembly, df_busco, df_annotation, assembly_summary_df):
    """Create a figure with summary tables and statistics"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Comprehensive Statistical Summary Tables', fontsize=16, fontweight='bold')
    
    # Assembly statistics table
    ax1 = axes[0, 0]
    if not df_assembly.empty and len(df_assembly.select_dtypes(include=[np.number]).columns) > 0:
        assembly_numeric = df_assembly.select_dtypes(include=[np.number])
        assembly_stats = assembly_numeric.describe()
        assembly_table_data = []
        for col in assembly_stats.columns:
            if col != 'Sample':
                assembly_table_data.append([
                    col.replace('_', ' '),
                    f"{assembly_stats.loc['mean', col]:.0f}",
                    f"{assembly_stats.loc['std', col]:.0f}",
                    f"{assembly_stats.loc['min', col]:.0f}",
                    f"{assembly_stats.loc['max', col]:.0f}"
                ])
        
        if assembly_table_data:
            table = ax1.table(cellText=assembly_table_data,
                            colLabels=['Metric', 'Mean', 'Std', 'Min', 'Max'],
                            cellLoc='center', loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
        else:
            ax1.text(0.5, 0.5, 'No assembly data available', ha='center', va='center', transform=ax1.transAxes)
    else:
        ax1.text(0.5, 0.5, 'No assembly data available', ha='center', va='center', transform=ax1.transAxes)
    
    ax1.set_title('Assembly Statistics', fontweight='bold')
    ax1.axis('off')
    
    # BUSCO statistics table
    ax2 = axes[0, 1]
    if not df_busco.empty and len(df_busco.select_dtypes(include=[np.number]).columns) > 0:
        busco_numeric = df_busco.select_dtypes(include=[np.number])
        busco_stats = busco_numeric.describe()
        busco_table_data = []
        for col in busco_stats.columns:
            if col != 'Sample':
                busco_table_data.append([
                    col.replace('_BUSCOs', '').replace('_', ' '),
                    f"{busco_stats.loc['mean', col]:.1f}%",
                    f"{busco_stats.loc['std', col]:.1f}%",
                    f"{busco_stats.loc['min', col]:.1f}%",
                    f"{busco_stats.loc['max', col]:.1f}%"
                ])
        
        if busco_table_data:
            table = ax2.table(cellText=busco_table_data,
                            colLabels=['BUSCO Category', 'Mean', 'Std', 'Min', 'Max'],
                            cellLoc='center', loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
        else:
            ax2.text(0.5, 0.5, 'No BUSCO data available', ha='center', va='center', transform=ax2.transAxes)
    else:
        ax2.text(0.5, 0.5, 'No BUSCO data available', ha='center', va='center', transform=ax2.transAxes)
    
    ax2.set_title('BUSCO Completeness Statistics', fontweight='bold')
    ax2.axis('off')
    
    # Annotation statistics table
    ax3 = axes[1, 0]
    if not df_annotation.empty and len(df_annotation.select_dtypes(include=[np.number]).columns) > 0:
        annotation_numeric = df_annotation.select_dtypes(include=[np.number])
        annotation_stats = annotation_numeric.describe()
        annotation_table_data = []
        for col in annotation_stats.columns:
            if col != 'Sample':
                annotation_table_data.append([
                    col,
                    f"{annotation_stats.loc['mean', col]:.0f}",
                    f"{annotation_stats.loc['std', col]:.0f}",
                    f"{annotation_stats.loc['min', col]:.0f}",
                    f"{annotation_stats.loc['max', col]:.0f}"
                ])
        
        if annotation_table_data:
            table = ax3.table(cellText=annotation_table_data,
                            colLabels=['Feature Type', 'Mean', 'Std', 'Min', 'Max'],
                            cellLoc='center', loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
        else:
            ax3.text(0.5, 0.5, 'No annotation data available', ha='center', va='center', transform=ax3.transAxes)
    else:
        ax3.text(0.5, 0.5, 'No annotation data available', ha='center', va='center', transform=ax3.transAxes)
    
    ax3.set_title('Annotation Feature Statistics', fontweight='bold')
    ax3.axis('off')
    
    # Pan-genome explanation
    ax4 = axes[1, 1]
    explanation_text = """
EXPLANATION OF KEY CONCEPTS:

PAN-GENOME STRUCTURE:
• Core Genes: Present in ≥99% of genomes
  - Essential genes shared by all/most isolates
  - Usually housekeeping and metabolic genes
  
• Shell Genes: Present in 15-99% of genomes  
  - Dispensable genes with intermediate frequency
  - Often adaptation-related functions
  
• Cloud Genes: Present in <15% of genomes
  - Rare/unique genes in few isolates
  - Often strain-specific or plasmid genes

GENE FREQUENCY DISTRIBUTION:
• X-axis: Number of genomes containing each gene family
• Y-axis: Number of gene families with that frequency
• Core Threshold (99% line): Genes present in almost all genomes
• Peak at right = many core genes
• Peak at left = many rare genes

ASSEMBLY QUALITY METRICS:
• N50: Length where 50% of assembly is in contigs ≥ this length
• Higher N50 = better assembly contiguity
• Fewer contigs = better assembly
• BUSCO completeness = gene content quality
    """
    
    ax4.text(0.05, 0.95, explanation_text, transform=ax4.transAxes,
            verticalalignment='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    ax4.set_title('Key Concepts & Interpretations', fontweight='bold')
    ax4.axis('off')
    
    plt.tight_layout()
    return fig

def main():
    """Main function for enhanced data visualization"""
    try:
        print("Creating comprehensive statistical analysis report...")
        
        # Create output directory
        output_dir = Path("results/enhanced_visualization")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Get sample list from config
        import yaml
        try:
            with open('config/config.yaml', 'r') as f:
                config = yaml.safe_load(f)
        except FileNotFoundError:
            print("ERROR: config/config.yaml not found")
            exit(1)
        
        if 'samples' not in config:
            print("ERROR: 'samples' not found in config")
            exit(1)
            
        try:
            samples_df = pd.read_csv(config['samples'], sep='\t')
            samples = samples_df['sample'].tolist()
        except Exception as e:
            print(f"ERROR: Could not read samples file: {e}")
            exit(1)
        
        print(f"Processing {len(samples)} samples...")
        
        # Load data from individual QUAST reports (fallback)
        print("Loading assembly statistics from QUAST reports...")
        assembly_data = []
        for sample in samples:
            quast_file = f"results/assembly_qc/{sample}/report.tsv"
            try:
                with open(quast_file, 'r') as f:
                    lines = f.readlines()
                    stats_dict = {'Sample': sample}
                    
                    for line in lines:
                        if '\t' in line:
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                metric = parts[0].strip()
                                value = parts[1].strip()
                                
                                if 'N50' in metric and 'N50' not in stats_dict:
                                    try:
                                        stats_dict['N50'] = int(value.replace(',', ''))
                                    except:
                                        pass
                                elif '# contigs' in metric:
                                    try:
                                        stats_dict['Contigs'] = int(value.replace(',', ''))
                                    except:
                                        pass
                                elif 'Total length' in metric:
                                    try:
                                        stats_dict['Total_length'] = int(value.replace(',', ''))
                                    except:
                                        pass
                                elif 'GC (%)' in metric:
                                    try:
                                        stats_dict['GC_content'] = float(value)
                                    except:
                                        pass
                    
                    assembly_data.append(stats_dict)
            except FileNotFoundError:
                print(f"Warning: QUAST report not found for {sample}")
                continue
        
        if not assembly_data:
            print("ERROR: No assembly data found")
            exit(1)
            
        df_assembly = pd.DataFrame(assembly_data)
        
        # Load assembly summary for boxplots
        print("Loading assembly summary statistics...")
        assembly_summary_df = load_assembly_summary_stats("results/assembly/assembly_summary.txt")
        print(f"Assembly summary DataFrame: {len(assembly_summary_df)} rows, columns: {list(assembly_summary_df.columns) if not assembly_summary_df.empty else 'empty'}")
        
        print("Loading BUSCO statistics...")
        df_busco = load_busco_stats(samples)
        print(f"BUSCO DataFrame: {len(df_busco)} rows, columns: {list(df_busco.columns) if not df_busco.empty else 'empty'}")
        
        print("Loading annotation statistics...")
        df_annotation = load_annotation_stats(samples)
        print(f"Annotation DataFrame: {len(df_annotation)} rows, columns: {list(df_annotation.columns) if not df_annotation.empty else 'empty'}")
        
        # Create comprehensive PDF report
        print("Creating comprehensive PDF report...")
        output_file = "results/enhanced_visualization/comprehensive_analysis_report.pdf"
        create_comprehensive_pdf_report(df_assembly, df_busco, df_annotation, assembly_summary_df, output_file)
        
        print("Comprehensive analysis report completed!")
        print(f"Report saved as: {output_file}")
        print("\nThis report includes:")
        print("  - Assembly quality distributions (boxplots, violin plots)")
        print("  - Assembly summary boxplots with statistics")
        print("  - Annotation feature analysis")
        print("  - Summary tables with descriptive statistics")
        print("  - Explanations of pan-genome concepts and metrics")
        
    except Exception as e:
        print(f"ERROR in visualization script: {e}")
        import traceback
        traceback.print_exc()
        exit(1)

if __name__ == "__main__":
    main()
