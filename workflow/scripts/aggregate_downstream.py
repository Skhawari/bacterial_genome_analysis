# scripts/aggregate_downstream.py

import pandas as pd
import glob
import yaml
from pathlib import Path

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
        print(f" {len(samples)} samples will be processed (excluded {len(blacklist)} blacklisted)\n")

    except Exception as e:
        print(f" Error reading sample or blacklist files: {e}")
        return

    # Output dir
    Path("results/downstream").mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter("results/downstream/summary.xlsx", engine="xlsxwriter") as writer:

        # --- Aggregation helper ---
        def aggregate_table(glob_path, sheet, filetype="tsv", infer_sample=False):
            files = glob.glob(glob_path)
            dfs = []
            for f in files:
                try:
                    df = pd.read_csv(f, sep="\t" if filetype == "tsv" else ",")
                    sample = Path(f).stem if not infer_sample else Path(f).parent.name
                    if sample not in samples:
                        continue
                    if infer_sample:
                        df["sample"] = sample
                    dfs.append(df)
                except Exception as e:
                    print(f"⚠️ Could not read {f}: {e}")
            if dfs:
                combined = pd.concat(dfs, ignore_index=True)
                combined.to_excel(writer, sheet_name=sheet, index=False)
                print(f" {sheet} → {len(combined)} entries")
            else:
                print(f"⚠️ No valid files found for {sheet}")

        # MLST
        if cfg["pipeline"].get("run_mlst", False):
            print(" Aggregating MLST...")
            aggregate_table("results/downstream/mlst/*.txt", "MLST", filetype="tsv")

        # AMR
        if cfg["pipeline"].get("run_amr", False):
            print(" Aggregating AMR...")
            aggregate_table("results/downstream/amr/*.tsv", "AMR")

        # Virulence
        if cfg["pipeline"].get("run_virulence", False):
            print(" Aggregating Virulence...")
            aggregate_table("results/downstream/virulence/*.tsv", "Virulence")

        # Plasmid (from blast hits)
        if cfg["pipeline"].get("run_plasmid", False):
            print(" Aggregating Plasmid BLAST...")
            aggregate_table("results/downstream/plasmid/*/blast_plasmid_hits.tsv", "Plasmid", infer_sample=True)

    print("\n Aggregation complete → results/downstream/summary.xlsx")

if __name__ == "__main__":
    main()