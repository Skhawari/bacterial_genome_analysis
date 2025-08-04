# Bacterial Genome Analysis Pipeline

This repository contains a Snakemake-based workflow for bacterial whole-genome sequencing (WGS) analysis. The pipeline supports both Illumina short reads and long reads (Oxford Nanopore or PacBio) and performs:

- Quality control (QC)
- De novo and hybrid genome assembly
- Assembly quality assessment (QUAST, BUSCO)
- Genome annotation (Prokka)
- Optional downstream analyses:
  - MLST typing
  - Antimicrobial resistance gene detection (RGI or ABRicate)
  - Virulence factor screening
  - Plasmid content analysis (BLAST + Mash)
- Comparative genomics and phylogeny (Roary + IQ-TREE)
- Final result aggregation and reporting

---

## 📁 Project Structure

```
.
├── results/                      # Output directory
├── config/
│   ├── config.yaml               # Main configuration file
│   ├── samples.tsv              # Sample sheet for short-read data
│   ├── samples_long_read.tsv    # Sample sheet for long-read data (optional)
│   └── blacklist.txt            # Optional: samples to exclude
├── workflow/                    
│   ├── envs/                         # Conda environment YAMLs
│   ├── rules/                        # Snakemake rule files
│   ├── scripts/                      # Custom helper scripts
│   ├── Snakefile                     # Entry point for Snakemake
└── README.md
```

---

## ⚙️ Requirements

- Linux/macOS
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Snakemake ≥ 7.0

---

## 🧪 Installation

```bash
# Clone this repository
git clone https://github.com/yourusername/bacterial-genome-pipeline.git
cd bacterial-genome-pipeline

# Create and activate base environment
conda create -n snakemake_env snakemake=7.32.4 python=3.11
conda activate snakemake_env
```

---

## 🚀 Running the Workflow

1. Place your raw FASTQ files in a folder of your choice.
2. Edit the sample sheet `config/samples.tsv` with sample names and paths.
3. Adjust parameters in `config/config.yaml` (e.g., enable/disable modules).
4. Run the workflow:

```bash
snakemake --cores 4
```

You can also generate the rule graph:

```bash
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

---

## ⚙️ Configuration

All pipeline options are defined in `config/config.yaml`:
- Paths to short- and long-read sample sheets
- Activation of optional modules (MLST, AMR, plasmid, etc.)
- Tool-specific parameters for SPAdes, Prokka, QUAST, BUSCO, Roary, IQ-TREE
- Blacklist and outgroup options

Lineage datasets and external databases must be downloaded manually or using provided Snakemake rules.

---

## 📊 Results

The output is structured in `results/` and includes:
- QC reports (FastQC, MultiQC)
- Assembly (FASTA) and quality assessments (QUAST, BUSCO)
- Annotation files (GFF, GBK, FASTA)
- Downstream results (TSV/Excel summaries)
- Phylogenetic trees (Newick)
- Visualizations and final reports

---

## 🧬 BUSCO & Plasmid Database Notes

- The BUSCO database (`bacteria_odb10`) is downloaded in a dedicated rule to avoid conflicts when running with multiple threads.
- The plasmid database (PLSDB) was manually downloaded to ensure compatibility with Mash and BLAST, due to issues with outdated tools not supporting the new format.

---

## 📄 Citation

If you use this pipeline in your research, please cite:

- [Snakemake](https://doi.org/10.1093/bioinformatics/bts480)
- [SPAdes](https://doi.org/10.1089/cmb.2012.0021)
- [QUAST](https://doi.org/10.1093/bioinformatics/btt086)
- [BUSCO](https://doi.org/10.1093/bioinformatics/btv351)
- [Prokka](https://doi.org/10.1093/bioinformatics/btu153)
- [Roary](https://doi.org/10.1093/bioinformatics/btv421)
- [IQ-TREE](https://doi.org/10.1093/molbev/msu300)

---

## 👩‍🔬 Authors

- Michael Riethmüller
- Sajjad Khawari 

Summer Term 2025 – Applied Sequence Analysis (Dr. Sandro Andreotti)
