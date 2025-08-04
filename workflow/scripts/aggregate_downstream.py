#!/usr/bin/env python3
import pandas as pd
import glob
from pathlib import Path
import os

def escape_excel_formulas(df):
    """Escape potential Excel formula characters for security"""
    for col in df.select_dtypes(include=['object']):
        if df[col].dtype == 'object':
            df[col] = df[col].astype(str).str.replace('\x00', '', regex=False)  # Remove null bytes
            df[col] = df[col].str.replace('^=', "'=", regex=True)   # Escape leading equals
            df[col] = df[col].str.replace('^\\+', "'+", regex=True) # Escape leading plus
            df[col] = df[col].str.replace('^-', "'-", regex=True)   # Escape leading minus
            df[col] = df[col].str.replace('^@', "'@", regex=True)   # Escape leading at
    return df

def read_abricate(module_name, pattern, samples):
    """Read and aggregate abricate output files (AMR/Virulence)"""
    dfs = []
    successful_files = 0
    
    for fp in glob.glob(pattern):
        try:
            sample = Path(fp).stem
            if sample not in samples:
                continue
            
            # Check if file exists and is not empty
            if not os.path.exists(fp) or os.path.getsize(fp) == 0:
                print(f" {sample}: Empty or missing file")
                continue
            
            # Read file with error handling
            df = pd.read_csv(fp, sep="\t", dtype=str, skiprows=1, header=None)
            if df.empty:
                print(f" {sample}: No data after header")
                continue
                
            # Add proper column names
            df.columns = [
                'FILE', 'SEQUENCE', 'START', 'END', 'STRAND', 'GENE', 'COVERAGE',
                'COVERAGE_MAP', 'GAPS', 'PERCENT_COVERAGE', 'PERCENT_IDENTITY',
                'DATABASE', 'ACCESSION', 'PRODUCT', 'RESISTANCE'
            ]
            df.insert(0, 'sample', sample)
            dfs.append(df)
            successful_files += 1
            print(f"   → {sample}: {len(df)} {module_name.lower()} hits")
            
        except Exception as e:
            print(f" Error reading {fp}: {e}")
            continue
    
    combined = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    print(f" {module_name} → {len(combined)} entries from {successful_files} files")
    return combined

def create_summary_report(downstream_cfg, samples):
    """Create summary report when no data is available"""
    summary_data = []
    
    for module in ['run_mlst', 'run_amr', 'run_virulence', 'run_plasmid']:
        if downstream_cfg.get(module, False):
            file_patterns = {
                'run_mlst': 'results/downstream/mlst/*.txt',
                'run_amr': 'results/downstream/amr/*.tsv',
                'run_virulence': 'results/downstream/virulence/*.tsv',
                'run_plasmid': 'results/downstream/plasmid/*/blast_plasmid_hits.tsv'
            }
            
            pattern = file_patterns[module]
            files_found = len(glob.glob(pattern))
            
            summary_data.append({
                'Module': module.replace('run_', '').upper(),
                'Enabled': True,
                'Files_Found': files_found,
                'Expected_Samples': len(samples),
                'Status': 'No data' if files_found == 0 else 'Files exist but empty or unreadable'
            })
    
    return pd.DataFrame(summary_data)

def main():
    print("Aggregating downstream analysis results...\n")
    
    # Load the Snakemake config directly
    cfg = snakemake.config
    downstream_cfg = cfg.get("downstream", {})

    # Collect sample list and apply complete_blacklist
    try:
        samples_df = pd.read_csv(cfg["samples"], sep="\t")
        samples = samples_df["sample"].tolist()
        
        # Apply blacklist
        bl_file = cfg.get("complete_blacklist", "")
        if bl_file and Path(bl_file).exists():
            with open(bl_file) as f:
                blacklist = [l.strip() for l in f if l.strip() and not l.startswith("#")]
            samples = [s for s in samples if s not in blacklist]
            print(f"Applied blacklist: {len(blacklist)} samples excluded")
        
        # Add outgroup if configured
        og_fasta = cfg.get("outgroup_fasta", "")
        if og_fasta:
            og_name = cfg.get("outgroup_name", "outgroup")
            if og_name not in samples:
                samples.append(og_name)
                print(f"Added outgroup: {og_name}")
                
        print(f"Processing {len(samples)} samples total\n")
        
    except Exception as e:
        print(f"Error reading sample configuration: {e}")
        return

    # Ensure output directory exists
    Path("results/downstream").mkdir(parents=True, exist_ok=True)

    # Collect all data with robust error handling
    all_data = {}
    modules_processed = 0

    # MLST Processing
    if downstream_cfg.get("run_mlst", False):
        print(" Aggregating MLST...")
        records = []
        successful_files = 0
        
        for fp in glob.glob("results/downstream/mlst/*.txt"):
            try:
                sample = Path(fp).stem
                if sample not in samples:
                    continue
                
                # Check file exists and is not empty
                if not os.path.exists(fp) or os.path.getsize(fp) == 0:
                    print(f" {sample}: Empty or missing MLST file")
                    continue
                
                # Read and parse MLST line
                lines = Path(fp).read_text().strip().splitlines()
                if not lines:
                    print(f" {sample}: Empty MLST file")
                    continue
                    
                line = lines[0].strip()
                parts = line.split("\t")
                
                if len(parts) >= 3:
                    records.append({
                        "sample": sample,
                        "file": parts[0],
                        "scheme": parts[1],
                        "sequence_type": parts[2],
                        "alleles": '\t'.join(parts[3:]) if len(parts) > 3 else ''
                    })
                    successful_files += 1
                    print(f" {sample}: {parts[1]} scheme, ST {parts[2]}")
                else:
                    print(f" {sample}: Invalid MLST format ({len(parts)} columns)")
                    
            except Exception as e:
                print(f"Error reading MLST file {fp}: {e}")
                continue
        
        if records:
            all_data["MLST"] = pd.DataFrame(records)
            modules_processed += 1
        print(f" MLST → {len(records)} entries from {successful_files} files\n")

    # AMR Processing
    if downstream_cfg.get("run_amr", False):
        print("Aggregating AMR...")
        amr_df = read_abricate("AMR", "results/downstream/amr/*.tsv", samples)
        if not amr_df.empty:
            all_data["AMR"] = amr_df
            modules_processed += 1
        print()

    # Virulence Processing
    if downstream_cfg.get("run_virulence", False):
        print("Aggregating Virulence...")
        vir_df = read_abricate("Virulence", "results/downstream/virulence/*.tsv", samples)
        if not vir_df.empty:
            all_data["Virulence"] = vir_df
            modules_processed += 1
        print()

    # Plasmid Processing
    if downstream_cfg.get("run_plasmid", False):
        print("Aggregating Plasmid BLAST...")
        dfs = []
        successful_files = 0
        
        for fp in glob.glob("results/downstream/plasmid/*/blast_plasmid_hits.tsv"):
            try:
                sample = Path(fp).parent.name
                if sample not in samples:
                    continue
                
                # Check file exists and is not empty
                if not os.path.exists(fp) or os.path.getsize(fp) == 0:
                    print(f" {sample}: Empty or missing plasmid file")
                    continue
                
                df = pd.read_csv(fp, sep="\t", header=None, dtype=str,
                               names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore'])
                if df.empty:
                    print(f" {sample}: No plasmid data")
                    continue
                    
                df.insert(0, 'sample', sample)
                dfs.append(df)
                successful_files += 1
                print(f" {sample}: {len(df)} plasmid hits")
                
            except Exception as e:
                print(f" Error reading plasmid file {fp}: {e}")
                continue
        
        combined = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
        if not combined.empty:
            all_data["Plasmid"] = combined
            modules_processed += 1
        print(f" Plasmid → {len(combined)} entries from {successful_files} files\n")

    # Write results with robust error handling and security
    excel_path = Path("results/downstream/summary.xlsx")
    
    if all_data:
        print(f" Writing results to Excel ({modules_processed} modules)...")
        
        # Remove existing file
        if excel_path.exists():
            excel_path.unlink()
        
        try:
            # Try to write Excel file with security measures
            with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
                for sheet_name, df in all_data.items():
                    # Ensure valid sheet name (Excel limits)
                    safe_name = sheet_name[:31].replace('[', '').replace(']', '').replace('*', '').replace('?', '').replace(':', '').replace('/', '_').replace('\\', '_')
                    
                    # Apply Excel security - escape formula characters
                    df_safe = escape_excel_formulas(df.copy())
                    
                    # Write to Excel
                    df_safe.to_excel(writer, sheet_name=safe_name, index=False)
                    print(f" Sheet '{safe_name}': {len(df_safe)} rows")
            
            print(f"\n Success! {excel_path} ({len(all_data)} sheets)")
            
        except Exception as e:
            print(f" Excel write failed: {e}")
            print(" Falling back to CSV files...")
            
            # Fallback mechanism - write CSV files
            for sheet_name, df in all_data.items():
                csv_file = Path(f"results/downstream/{sheet_name.lower()}_summary.csv")
                df_safe = escape_excel_formulas(df.copy())
                df_safe.to_csv(csv_file, index=False)
                print(f" Wrote CSV: {csv_file}")
            
            print(f"\n Fallback complete → {len(all_data)} CSV files in results/downstream/")
    
    else:
        print("No downstream data found to aggregate.")
        print("Creating summary report...")
        
        # Create summary report even when no data
        summary_df = create_summary_report(downstream_cfg, samples)
        
        if not summary_df.empty:
            try:
                with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
                    summary_df.to_excel(writer, sheet_name="Summary", index=False)
                print(f" Created summary report → {excel_path}")
            except Exception as e:
                print(f" Could not write summary: {e}")
        else:
            print(" No downstream modules enabled in config")

    print("\n Aggregation completed!")

if __name__ == "__main__":
    main()
