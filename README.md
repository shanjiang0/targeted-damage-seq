# Targeted Damage-seq Processing Pipeline

Processing pipeline for **Targeted Damage-seq** paired-end sequencing data, from raw FASTQ files to damage-site BED files, coverage bedGraphs, and nucleotide/k-mer composition at damage sites.

This code accompanies the manuscript *Targeted DNA damage sequencing and multi-omic analysis reveal repair kinetics and genomic context shaping UV-induced mutagenesis* (under review).

> **Note:** This repository contains only the Damage-seq **processing pipeline**. Downstream analysis and figure-reproduction code for the paper is available at [UV-tDamage-seq-mutagenesis](https://github.com/shanjiang0/UV-tDamage-seq-mutagenesis) (Zenodo DOI to be added upon publication).

---

## Repository Structure

```
.
├── README.md
├── LICENSE
├── .gitignore
└── pipeline/
    ├── build_project.sh                     # Sets up folders and generates all job scripts
    ├── Targeted_Damageseq_PE_pipeline.py    # Main pipeline (paired-end)
    ├── Example/
    │   └── 00.rawdata/
    │       └── targeted_Damageseq_sample.list   # Example sample list
    └── utils/                               # Helper scripts called by the pipeline
        ├── format_unify_fastq_check.py
        ├── filter_read1_PE.py
        ├── Dedup_PE.sh
        ├── fix_bed.py
        ├── fa2nt_abundance.py
        └── kmer_abundance.py
```

---

## Pipeline Steps

| Step | Description | Tool |
|---|---|---|
| 01 | Merge / link FASTQ files per sample | `cat` / `ln` |
| 02 | Adapter trimming | cutadapt |
| 03 | FASTQ formatting and quality filtering (Q ≥ 20) | Python / Biopython |
| 04 | Deduplication | `Dedup_PE.sh` |
| 05 | Alignment | Bowtie1 |
| 06 | SAM → BAM → BED (properly paired, Read 1 only) | SAMtools, bedtools |
| 07 | Replicate merging | `cat` |
| 08 | Slop, range fixing, bedGraph generation | bedtools |
| 09 | FASTA extraction, nucleotide and k-mer abundance | bedtools, Python |

---

## Requirements

The pipeline is designed for a **Linux HPC cluster with SLURM** and environment modules (Lmod). Generated job scripts load these modules:

| Software | Tested version |
|---|---|
| cutadapt | 4.9 |
| Biopython | 1.84 |
| Bowtie | 1.3.1 |
| SAMtools | 1.21 |
| bedtools | 2.31.0 |
| Python | 3.12 |

If your system does not use modules, install these tools yourself and remove the `ml ...` lines from the generated scripts. Module names may also differ on your cluster.

### Reference files

- **Bowtie1 index** built from your reference genome (e.g. hg19 or hg38)
- **Reference FASTA**, the same genome used to build the index
- **Chromosome sizes file** (UCSC format) for the same genome

---

## Usage

### 1. Prepare the sample list

The sample list is a **tab-separated** file with 3 columns and **no header line**:

| Column | Description | Example |
|---|---|---|
| 1. rep_name | Replicate group name. Samples with the same rep_name are merged in step 07. | `1` |
| 2. sampleID | Unique sample ID | `CPD1` |
| 3. fq1,fq2 | Paths to Read 1 and Read 2 FASTQ files, separated by a comma (no space) | `/path/R1.fastq.gz,/path/R2.fastq.gz` |

See the example: [`pipeline/Example/00.rawdata/targeted_Damageseq_sample.list`](pipeline/Example/00.rawdata/targeted_Damageseq_sample.list)

The file must be named `targeted_Damageseq_sample.list` and placed in `00.rawdata/` inside your working directory. Columns must be separated by real tabs, not spaces.

### 2. Configure and run the setup script

```bash
mkdir my_project && cd my_project
mkdir -p 00.rawdata
cp /path/to/your/targeted_Damageseq_sample.list 00.rawdata/
cp /path/to/target-damage-seq/pipeline/build_project.sh .
```

Edit the **Configuration** section of `build_project.sh`:

```bash
treatment="UV"        # UV, Cisplatin, or untreated
bt_ref="/path/to/bowtie1_index/hg19"
ref_fa="/path/to/reference/hg19.fa"
chrom_sizes_slop="/path/to/reference/hg19.chrom.sizes"
pipeline_script="/path/to/target-damage-seq/pipeline/Targeted_Damageseq_PE_pipeline.py"
```

Then run:

```bash
bash build_project.sh
```

This creates the folder structure and writes one job script per sample and per replicate group.

### 3. Submit the jobs

Run the preprocessing jobs first, and wait for them to finish:

```bash
cd 01.scripts/01.preprocessing && sbatch sub_array.sh
```

Then run replicate merging:

```bash
cd ../02.replicates_merging && sbatch sub_array.sh
```

Adjust the `#SBATCH` settings (partition, memory, time) in the generated `sub_array.sh` files for your cluster.

---

## Output

```
my_project/
├── 00.rawdata/                      # Sample list
├── 01.scripts/                      # Generated job scripts
│   ├── 01.preprocessing/
│   └── 02.replicates_merging/
└── 02.results/
    ├── 01.preprocessing/<sampleID>/ # Trimmed/filtered FASTQ, BAM, <sampleID>.sorted.bed
    └── 02.replicates_merging/<rep_name>/
```

Key output files per replicate group (`02.results/02.replicates_merging/<rep_name>/`):

| File | Description |
|---|---|
| `<rep>.rep.merged.bed` | Merged damage-site reads from all samples in the group |
| `<rep>.combined2_sorted.bedgraph` | Coverage of 2-bp damage sites |
| `<rep>_nt_abun.txt`, `<rep>_kmer_abun.txt` | Nucleotide and dinucleotide composition, 10-bp window |
| `<rep>_nt_abun2.txt`, `<rep>_kmer_abun2.txt` | Nucleotide and dinucleotide composition at the 2-bp damage site |

---

## Citation

The associated manuscript is currently **under review**. Until it is published, please cite this repository:

> Jiang, S. (2026). *TARGET Damage-seq Processing Pipeline* (Version 1.0.0) [Computer software]. https://github.com/shanjiang0/target-damage-seq

You can also click **"Cite this repository"** on the right side of this GitHub page for APA and BibTeX formats.

The manuscript citation will be added here upon publication:

> Jiang, S., *et al.* Targeted DNA damage sequencing and multi-omic analysis reveal repair kinetics and genomic context shaping UV-induced mutagenesis. *Manuscript under review.*

---

## Author and contact

Code written by **Shan Jiang** (GitHub: [@shanjiang0](https://github.com/shanjiang0)), Li Lab, University of Georgia.

Corresponding author: **Dr. Wentao Li** (Wentao.Li@uga.edu), University of Georgia.

For questions about the code, please open a [GitHub Issue](../../issues) or contact the corresponding author.

---

## License

This project is licensed under the [MIT License](LICENSE).
