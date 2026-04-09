# Targeted DNA Damage Sequencing and Integrative Analysis Identify Genomic Determinants of Melanoma Mutagenesis

> **Note:** This repository contains code associated with the manuscript currently under review. The primary dataset will be made publicly available upon publication.

---

## Overview

This repository provides the full analysis pipeline for **TARGET Damage-seq**, including data processing, statistical analysis, machine learning models, and figure reproduction scripts used in the paper. The pipeline integrates our newly generated Damage-seq data with published multi-omics datasets to identify genomic determinants of melanoma mutagenesis.

---

## Repository Structure

```
.
├── README.md
├── envs/                        # Conda/environment configuration files
│   ├── environment_python.yml
│   └── environment_r.yml
├── data/
│   ├── raw/                     # Raw input data (not included; see Data Availability)
│   └── processed/               # Intermediate processed files
├── pipeline/
│   └── damage_seq_pipeline.sh   # Main Bash pipeline for Damage-seq processing
├── analysis/
│   ├── preprocessing.py         # Data cleaning and feature engineering
│   ├── integrative_analysis.py  # Multi-omics integration and ML models
│   └── statistics.py            # Statistical tests and models
├── R/
│   ├── visualization.R          # Main figure plotting scripts
│   ├── genomic_analysis.R       # Genomic context and mutation signature analysis
│   └── differential_analysis.R # Differential damage analysis
├── notebooks/                   # (Optional) Exploratory Jupyter / R notebooks
├── results/
│   └── figures/                 # Output figures
└── scripts/                     # Miscellaneous utility scripts
```

> **Note:** Adjust folder names to match your actual structure before pushing.

---

## Data Availability

| Dataset | Source | Accession |
|---|---|---|
| TARGET Damage-seq (this study) | GEO | Pending release upon publication |
| Published omics datasets | GEO / public repositories | See paper Methods section |

Once the paper is published, the primary dataset will be released. In the meantime, contact the corresponding author for data access inquiries.

---

## Requirements

### System
- Linux / macOS (Bash pipeline)
- Python ≥ 3.8
- R ≥ 4.1

### Python Dependencies
Install via conda or pip:
```bash
conda env create -f envs/environment_python.yml
conda activate damage-seq-py
```
Or with pip:
```bash
pip install -r requirements.txt
```

Key Python packages: `numpy`, `pandas`, `scikit-learn`, `scipy`, `matplotlib`, `seaborn`

### R Dependencies
```r
# In R console
install.packages(c("ggplot2", "dplyr", "tidyr", "BiocManager"))
BiocManager::install(c("GenomicRanges", "BSgenome", "SomaticSignatures"))
# Add other Bioconductor packages as needed
```

---

## Usage

### 1. Damage-seq Processing Pipeline (Bash)
```bash
bash pipeline/damage_seq_pipeline.sh \
  --input data/raw/sample.fastq.gz \
  --genome hg38 \
  --output data/processed/
```

### 2. Integrative Analysis (Python)
```bash
python analysis/integrative_analysis.py \
  --damage data/processed/damage_sites.bed \
  --omics data/raw/omics/ \
  --output results/
```

### 3. Figure Generation (R)
```r
Rscript R/visualization.R
```

> See individual script headers for full argument documentation.

---

## Reproducing Paper Results

To reproduce the main figures and analyses from the paper:

1. Download data (see **Data Availability** above)
2. Place raw files under `data/raw/`
3. Run the Bash pipeline to generate processed files
4. Run Python analysis scripts
5. Run R visualization scripts to generate figures under `results/figures/`

A step-by-step walkthrough is provided in `notebooks/` (coming soon).

---

## Citation

If you use this code or data in your work, please cite:

> [Authors]. *Targeted DNA Damage Sequencing and Integrative Analysis Identify Genomic Determinants of Melanoma Mutagenesis.* [Journal], [Year]. [DOI — to be added upon publication]

---

## License

This project is licensed under the [MIT License](LICENSE) — feel free to use and adapt with attribution.

---

## Contact

For questions about the code or data, please open a [GitHub Issue](../../issues) or contact the corresponding author at **[your.email@institution.edu]**.
