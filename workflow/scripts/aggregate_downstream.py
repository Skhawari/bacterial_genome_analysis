# scripts/aggregate_downstream.py

import pandas as pd
import glob
import yaml
from pathlib import Path
import sys
import os

def main():
    print("🔎 Aggregating downstream analysis results...\n")

    # Load config
    try:
        with open("config/config.yaml", 'r') as f:
            cfg = yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
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
        if cfg["annotation"].get("outgroup_fasta", ""):
            outgroup_name = cfg["annotation"].get("outgroup_name", "outgroup")
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
                    # Read with header starting with # - skip comment parameter
                    df = pd.read_csv(f, sep="\t", dtype=str)
                    # Remove # from column names if present
                    df.columns = df.columns.str.replace('^#', '', regex=True)
                    if not df.empty:
                        df.insert(0, 'sample', sample)
                        amr_data.append(df)
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
                    # Read with header starting with # - skip comment parameter
                    df = pd.read_csv(f, sep="\t", dtype=str)
                    # Remove # from column names if present
                    df.columns = df.columns.str.replace('^#', '', regex=True)
                    if not df.empty:
                        df.insert(0, 'sample', sample)
                        vir_data.append(df)
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
            except Exception as e:
                print(f"⚠️ Could not read plasmid file {f}: {e}")
                
        if plasmid_data:
            all_data['Plasmid'] = pd.concat(plasmid_data, ignore_index=True)
            print(f" Plasmid → {len(all_data['Plasmid'])} entries")

    # Write Excel file only if we have data
    if all_data:
        try:
            # Remove any existing corrupted file
            excel_path = "results/downstream/summary.xlsx"
            if os.path.exists(excel_path):
                os.remove(excel_path)
                
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
            
        except Exception as e:
            print(f"❌ Error writing Excel file: {e}")
            # Fallback: write as CSV files
            for sheet_name, df in all_data.items():
                csv_file = f"results/downstream/{sheet_name.lower()}_summary.csv"
                df.to_csv(csv_file, index=False)
                print(f"📄 Wrote fallback CSV: {csv_file}")
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