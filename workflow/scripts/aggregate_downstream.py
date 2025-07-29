# scripts/aggregate_downstream.py

import pandas as pd
import glob
import yaml
from pathlib import Path
import sys
import os

def main():
    print("🔎 Aggregating downstream analysis results...\n")

    # Load config - try to detect which config file is being used
    config_files = ["config/config_set2.yaml", "config/config.yaml"]
    cfg = None
    
    for config_file in config_files:
        try:
            if Path(config_file).exists():
                with open(config_file, 'r') as f:
                    cfg = yaml.safe_load(f)
                print(f"📄 Using config: {config_file}")
                break
        except Exception as e:
            print(f"⚠️ Could not load {config_file}: {e}")
            continue
    
    if cfg is None:
        print("❌ No valid config file found")
        return

    # Get samples & blacklist
    try:
        samples_df = pd.read_csv(cfg["samples"], sep="\t")
        all_samples = samples_df["sample"].tolist()

        blacklist = []
        if cfg.get("complete_blacklist") and Path(cfg["complete_blacklist"]).exists():
            with open(cfg["complete_blacklist"], 'r') as f:
                blacklist = [line.strip() for line in f if line.strip() and not line.startswith('#')]

        samples = [s for s in all_samples if s not in blacklist]
        
        # Add outgroup if configured
        if cfg.get("outgroup_fasta", ""):
            outgroup_name = cfg.get("outgroup_name", "outgroup")
            if outgroup_name not in samples:
                samples.append(outgroup_name)
                
        print(f" {len(samples)} samples will be processed (excluded {len(blacklist)} blacklisted)\n")

    except Exception as e:
        print(f" Error reading sample or blacklist files: {e}")
        return

    # Output dir
    Path("results/downstream").mkdir(parents=True, exist_ok=True)

    # Collect all data first to ensure we have something to write
    all_data = {}
    
    # MLST Processing (custom format)
    if cfg["pipeline"].get("run_mlst", False):
        print(" Aggregating MLST...")
        mlst_files = glob.glob("results/downstream/mlst/*.txt")
        mlst_data = []
        
        for f in mlst_files:
            try:
                sample = Path(f).stem
                if sample not in samples:
                    continue
                    
                with open(f, 'r') as file:
                    line = file.readline().strip()
                    if line:
                        # MLST format: file\tscheme\tST\tgene1\tgene2...
                        parts = line.split('\t')
                        if len(parts) >= 3:
                            mlst_data.append({
                                'sample': sample,
                                'file': parts[0],
                                'scheme': parts[1],
                                'sequence_type': parts[2],
                                'alleles': '\t'.join(parts[3:]) if len(parts) > 3 else ''
                            })
                            print(f"   → {sample}: {parts[1]} scheme, ST {parts[2]}")
            except Exception as e:
                print(f"⚠️ Could not read MLST file {f}: {e}")
                
        if mlst_data:
            all_data['MLST'] = pd.DataFrame(mlst_data)
            print(f" MLST → {len(mlst_data)} entries")

    # AMR Processing
    if cfg["pipeline"].get("run_amr", False):
        print(" Aggregating AMR...")
        amr_files = glob.glob("results/downstream/amr/*.tsv")
        amr_data = []
        
        for f in amr_files:
            try:
                sample = Path(f).stem
                if sample not in samples:
                    continue
                    
                if os.path.getsize(f) > 0:  # Check if file is not empty
                    # AMR files have problematic headers - read without header and add manually
                    try:
                        # Skip the header line and read data
                        df = pd.read_csv(f, sep="\t", dtype=str, skiprows=1, header=None)
                        if not df.empty:
                            # Add column names manually based on AMR output format
                            df.columns = ['FILE', 'SEQUENCE', 'START', 'END', 'STRAND', 'GENE', 'COVERAGE', 'COVERAGE_MAP', 'GAPS', 'PERCENT_COVERAGE', 'PERCENT_IDENTITY', 'DATABASE', 'ACCESSION', 'PRODUCT', 'RESISTANCE']
                            df.insert(0, 'sample', sample)
                            amr_data.append(df)
                            print(f"   → {sample}: {len(df)} AMR hits")
                    except Exception as parse_error:
                        print(f"⚠️ Parse error for AMR file {f}: {parse_error}")
            except Exception as e:
                print(f"⚠️ Could not read AMR file {f}: {e}")
                
        if amr_data:
            all_data['AMR'] = pd.concat(amr_data, ignore_index=True)
            print(f" AMR → {len(all_data['AMR'])} entries")

    # Virulence Processing
    if cfg["pipeline"].get("run_virulence", False):
        print(" Aggregating Virulence...")
        vir_files = glob.glob("results/downstream/virulence/*.tsv")
        vir_data = []
        
        for f in vir_files:
            try:
                sample = Path(f).stem
                if sample not in samples:
                    continue
                    
                if os.path.getsize(f) > 0:  # Check if file is not empty
                    # Virulence files have same header format as AMR - read without header
                    try:
                        # Skip the header line and read data
                        df = pd.read_csv(f, sep="\t", dtype=str, skiprows=1, header=None)
                        if not df.empty:
                            # Add column names manually based on virulence output format
                            df.columns = ['FILE', 'SEQUENCE', 'START', 'END', 'STRAND', 'GENE', 'COVERAGE', 'COVERAGE_MAP', 'GAPS', 'PERCENT_COVERAGE', 'PERCENT_IDENTITY', 'DATABASE', 'ACCESSION', 'PRODUCT', 'RESISTANCE']
                            df.insert(0, 'sample', sample)
                            vir_data.append(df)
                            print(f"   → {sample}: {len(df)} virulence factors")
                    except Exception as parse_error:
                        print(f"⚠️ Parse error for virulence file {f}: {parse_error}")
            except Exception as e:
                print(f"⚠️ Could not read virulence file {f}: {e}")
                
        if vir_data:
            all_data['Virulence'] = pd.concat(vir_data, ignore_index=True)
            print(f" Virulence → {len(all_data['Virulence'])} entries")

    # Plasmid Processing
    if cfg["pipeline"].get("run_plasmid", False):
        print(" Aggregating Plasmid BLAST...")
        plasmid_files = glob.glob("results/downstream/plasmid/*/blast_plasmid_hits.tsv")
        plasmid_data = []
        
        for f in plasmid_files:
            try:
                sample = Path(f).parent.name
                if sample not in samples:
                    continue
                    
                if os.path.getsize(f) > 0:  # Check if file is not empty
                    # BLAST format: qseqid sseqid pident length evalue bitscore
                    df = pd.read_csv(f, sep="\t", dtype=str, header=None, 
                                   names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore'])
                    if not df.empty:
                        df.insert(0, 'sample', sample)
                        plasmid_data.append(df)
                        print(f"   → {sample}: {len(df)} plasmid hits")
            except Exception as e:
                print(f"⚠️ Could not read plasmid file {f}: {e}")
                
        if plasmid_data:
            all_data['Plasmid'] = pd.concat(plasmid_data, ignore_index=True)
            print(f" Plasmid → {len(all_data['Plasmid'])} entries")

    # Write Excel file only if we have data
    if all_data:
        try:
            # Try to write Excel file
            excel_path = "results/downstream/summary.xlsx"
            if os.path.exists(excel_path):
                os.remove(excel_path)
                
            try:
                with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
                    for sheet_name, df in all_data.items():
                        # Ensure valid sheet name (Excel limit 31 chars, no special chars)
                        safe_name = sheet_name[:31].replace('[', '').replace(']', '').replace('*', '').replace('?', '').replace(':', '').replace('/', '_').replace('\\', '_')
                        
                        # Clean data for Excel compatibility - CRITICAL: Handle formula characters
                        for col in df.columns:
                            if df[col].dtype == 'object':
                                # Replace characters that Excel interprets as formulas
                                df[col] = df[col].astype(str).str.replace('\x00', '', regex=False)  # Remove null bytes
                                df[col] = df[col].str.replace('^=', "'=", regex=True)  # Escape leading equals
                                df[col] = df[col].str.replace('^\\+', "'+", regex=True)  # Escape leading plus
                                df[col] = df[col].str.replace('^-', "'-", regex=True)  # Escape leading minus
                                df[col] = df[col].str.replace('^@', "'@", regex=True)  # Escape leading at
                        
                        df.to_excel(writer, sheet_name=safe_name, index=False)
                        
                print(f"\n✅ Aggregation complete → results/downstream/summary.xlsx ({len(all_data)} sheets)")
                
            except ImportError:
                print("⚠️ openpyxl not available, writing CSV files instead...")
                raise Exception("openpyxl not found")
                
        except Exception as e:
            print(f"⚠️ Excel write failed: {e}")
            # Fallback: write as CSV files
            for sheet_name, df in all_data.items():
                csv_file = f"results/downstream/{sheet_name.lower()}_summary.csv"
                df.to_csv(csv_file, index=False)
                print(f"📄 Wrote CSV: {csv_file}")
            print(f"\n✅ Aggregation complete → {len(all_data)} CSV files in results/downstream/")
    else:
        print("⚠️ No data found to aggregate. Creating summary report.")
        # Create a summary of what was attempted
        summary_data = []
        for module in ['run_mlst', 'run_amr', 'run_virulence', 'run_plasmid']:
            if cfg["pipeline"].get(module, False):
                file_pattern = {
                    'run_mlst': 'results/downstream/mlst/*.txt',
                    'run_amr': 'results/downstream/amr/*.tsv',
                    'run_virulence': 'results/downstream/virulence/*.tsv',
                    'run_plasmid': 'results/downstream/plasmid/*/blast_plasmid_hits.tsv'
                }[module]
                
                files_found = len(glob.glob(file_pattern))
                summary_data.append({
                    'Module': module.replace('run_', '').upper(),
                    'Enabled': True,
                    'Files_Found': files_found,
                    'Status': 'No data' if files_found == 0 else 'Files exist but empty or unreadable'
                })
        
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            with pd.ExcelWriter("results/downstream/summary.xlsx", engine="openpyxl") as writer:
                summary_df.to_excel(writer, sheet_name="Summary", index=False)
            print("📄 Created summary report")
        else:
            print("⚠️ No downstream modules enabled")

if __name__ == "__main__":
    main()