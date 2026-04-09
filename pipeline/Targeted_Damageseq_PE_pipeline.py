#!/bin/python
"""
TARGET Damage-seq Paired-End Processing Pipeline
-------------------------------------------------
Usage:
    python Damageseq_PE_pipeline.py <sample_list> <work_dir> <treatment> <bt_ref> <ref_fa> <chrom_sizes_slop>

Arguments:
    sample_list      : TSV file with 3 columns: rep_name  sampleID  fq1,fq2
    work_dir         : Working directory for outputs
    treatment        : Treatment type (UV, Cisplatin, untreated)
    bt_ref           : Path to Bowtie1 index prefix (e.g. /path/to/bowtie1/hg38)
    ref_fa           : Path to reference FASTA file (e.g. /path/to/hg38.fa)
    chrom_sizes_slop : Path to chromosome sizes file for bedtools slop

Pipeline Steps:
    01. Merge/link FASTQ files per sample
    02. Adapter trimming (cutadapt)
    03. FASTQ formatting and QC filtering (Q>=20)
    04. Deduplication
    05. Alignment (Bowtie1)
    06. SAM -> BAM -> BED conversion
    07. Replicate merging
    08. Slop, range fixing, bedgraph generation
    09. FASTA extraction and nucleotide/kmer abundance

Dependencies:
    cutadapt/4.9, Biopython/1.84, Bowtie/1.3.1, SAMtools/1.21, pybedtools/0.10.0

Notes:
    - Bowtie1 index must be built from the same reference FASTA
    - Chromosome sizes file should match the reference genome (UCSC format)
    - utils/ scripts must be in the same directory as this pipeline script
"""

import sys
import os
from collections import defaultdict

# ── Configuration ──────────────────────────────────────────────────────────────

if len(sys.argv) != 7:
    sys.exit(
        "Usage: python Damageseq_PE_pipeline.py "
        "<sample_list> <work_dir> <treatment> <bt_ref> <ref_fa> <chrom_sizes_slop>"
    )

METHOD           = "Damageseq"
script_path      = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils/")
sample_list      = sys.argv[1]
work_dir         = sys.argv[2]
treatment        = sys.argv[3]
bt_ref           = sys.argv[4]
ref              = sys.argv[5]
chrom_sizes_slop = sys.argv[6]

# ── Bowtie1 parameters by treatment ───────────────────────────────────────────

bowtie_param = {
    "UV":        "-q --phred33-quals --nomaqround -m 4 --tryhard --allow-contain -p 4 --seed=123 -X 1000",
    "Cisplatin": "-q --phred33-quals --nomaqround -m 4 --tryhard --allow-contain -p 4 --seed=123 -X 1000",
    "untreated": "-q --phred33-quals --nomaqround -m 4 --tryhard --allow-contain -p 4 --seed=123 -X 1000",
}
bowtie_default = "-x -q --nomaqround -m 4 -v 3 --tryhard --strata --best -p 4 --seed=123"

# ── Step 1: Parse sample list ──────────────────────────────────────────────────

sample_dict  = dict()
merged_dict  = defaultdict(list)   # rep_name -> [sampleIDs]

with open(sample_list) as F1:
    for line in F1:
        rep_name, sampleID, FQs = line.strip().split("\t")
        fq1, fq2 = FQs.split(",")
        if sampleID not in sample_dict:
            sample_dict[sampleID] = {"rep_name": rep_name, "fq1": [fq1], "fq2": [fq2]}
            merged_dict[rep_name].append(sampleID)
        else:
            if rep_name != sample_dict[sampleID]["rep_name"]:
                sys.exit(f"ERROR: {sampleID} has conflicting replicate names.")
            sample_dict[sampleID]["fq1"].append(fq1)
            sample_dict[sampleID]["fq2"].append(fq2)

# ── Step 2: Write per-sample preprocessing scripts ────────────────────────────

script_list = f"{work_dir}/01.scripts/01.preprocessing/script.list"
script_num  = 0

with open(script_list, "w") as S1:
    for sampleID, val in sample_dict.items():
        rep_name   = val["rep_name"]
        fq1_list   = val["fq1"]
        fq2_list   = val["fq2"]
        sample_dir = f"{work_dir}/02.results/01.preprocessing/{sampleID}"
        run_sh     = f"{work_dir}/01.scripts/01.preprocessing/{sampleID}_run.sh"
        script_num += 1
        S1.write(f"{run_sh}\n")

        with open(run_sh, "w") as W1:

            # Create sample directory
            W1.write(f"mkdir -p {sample_dir}\n\n")

            # Merge or link FASTQ R1/R2
            out_fq1 = f"{sample_dir}/{sampleID}.merged_R1.fq.gz"
            out_fq2 = f"{sample_dir}/{sampleID}.merged_R2.fq.gz"
            W1.write("cat %s > %s\n" % (" ".join(fq1_list), out_fq1) if len(fq1_list) > 1
                     else "ln -s %s %s\n" % (fq1_list[0], out_fq1))
            W1.write("cat %s > %s\n" % (" ".join(fq2_list), out_fq2) if len(fq2_list) > 1
                     else "ln -s %s %s\n" % (fq2_list[0], out_fq2))

            # Adapter trimming
            cutadapt_param = (
                "-g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT "
                "--pair-filter=any --discard-trimmed "
                "--minimum-length=10 --no-indels -e 0.067 --overlap=6"
            )
            cut_fq1 = f"{sample_dir}/{sampleID}_merged_cutadapt_R1.fq.gz"
            cut_fq2 = f"{sample_dir}/{sampleID}_merged_cutadapt_R2.fq.gz"
            W1.write(f"\nml cutadapt/4.9-GCCcore-13.3.0\n")
            W1.write(f"cutadapt {cutadapt_param} -o {cut_fq1} -p {cut_fq2} {out_fq1} {out_fq2}\n")

            # FASTQ formatting
            format_fq1 = f"{sample_dir}/{sampleID}_cutadapt_format_R1.fq.gz"
            format_fq2 = f"{sample_dir}/{sampleID}_cutadapt_format_R2.fq.gz"
            W1.write(f"\npython3 {script_path}format_unify_fastq_check.py {cut_fq1} {cut_fq2} {format_fq1} {format_fq2}\n")

            # Quality filtering (Q >= 20)
            filtered_fq1 = f"{sample_dir}/{sampleID}_cutadapt_format_filtered_R1.fq.gz"
            filtered_fq2 = f"{sample_dir}/{sampleID}_cutadapt_format_filtered_R2.fq.gz"
            W1.write(f"\nml Biopython/1.84-gfbf-2024a\n")
            W1.write(f"python3 {script_path}filter_read1_PE.py {format_fq1} {format_fq2} {filtered_fq1} {filtered_fq2} 20\n")

            # Deduplication
            dedup_fq1 = f"{sample_dir}/{sampleID}_cutadapt_format_filtered_dedup_R1.fq"
            dedup_fq2 = f"{sample_dir}/{sampleID}_cutadapt_format_filtered_dedup_R2.fq"
            W1.write(f"\nsh {script_path}Dedup_PE.sh {filtered_fq1} {filtered_fq2} {dedup_fq1} {dedup_fq2}\n")

            # Bowtie1 alignment
            bt_param  = bowtie_param.get(treatment, bowtie_default)
            mapped_sam = f"{sample_dir}/{sampleID}_cutadapt_format_filtered_mapped.sam"
            W1.write(f"\nml Bowtie/1.3.1-GCC-13.3.0\n")
            W1.write(f"bowtie {bt_param} -x {bt_ref} -1 {dedup_fq1}.gz -2 {dedup_fq2}.gz -S {mapped_sam}\n")

            # SAM -> BAM -> read1 -> sorted
            mapped_bam  = f"{sample_dir}/{sampleID}_cutadapt_format_filtered_mapped.bam"
            filtered_bam= f"{sample_dir}/{sampleID}_cutadapt_format_filtered_mapped_properly_paired.bam"
            read1_bam   = f"{sample_dir}/{sampleID}_cutadapt_format_filtered_mapped_properly_paired_read1.bam"
            sorted_bam  = f"{sample_dir}/{sampleID}_cutadapt_format_filtered_mapped_properly_paired_sorted.bam"
            W1.write(f"\nml SAMtools/1.21-GCC-13.3.0\n")
            W1.write(f"samtools view -b -S -o {mapped_bam} {mapped_sam}\n")
            W1.write(f"samtools view -h -f 0x2 -b {mapped_bam} > {filtered_bam}\n")
            W1.write(f"samtools view -h -f 64  -b {filtered_bam} > {read1_bam}\n")
            W1.write(f"samtools sort {read1_bam} -o {sorted_bam}\n")

            # BAM -> BED
            sorted_bed = f"{sample_dir}/{sampleID}.sorted.bed"
            W1.write(f"\nml pybedtools/0.10.0-foss-2023b\n")
            W1.write(f"bedtools bamtobed -i {sorted_bam} | sort -u -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n > {sorted_bed}\n")

# SLURM array submission script for Step 1
with open(f"{work_dir}/01.scripts/01.preprocessing/sub_array.sh", "w") as O1:
    O1.write('\n'.join([
        '#!/bin/bash',
        f'#SBATCH --job-name={METHOD}_{treatment}_preprocessing_array',
        '#SBATCH --partition=batch',
        '#SBATCH --ntasks=1',
        '#SBATCH --cpus-per-task=64',
        '#SBATCH --time=120:00:00',
        '#SBATCH --mem=64G',
        '#SBATCH --output=%x_%j.out',
        '#SBATCH --mail-type=END,FAIL',
        f'#SBATCH --array=1-{script_num}',
        'file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" script.list)',
        'sh $file',
    ]) + "\n")

# ── Step 3: Write per-replicate merging scripts ────────────────────────────────

script_list2 = f"{work_dir}/01.scripts/02.replicates_merging/script.list"

with open(script_list2, "w") as S2:
    for rep_name, sample_ids in merged_dict.items():
        rep_dir    = f"{work_dir}/02.results/02.replicates_merging/{rep_name}"
        rep_run_sh = f"{work_dir}/01.scripts/02.replicates_merging/{rep_name}_run.sh"
        S2.write(f"{rep_run_sh}\n")

        with open(rep_run_sh, "w") as W2:
            W2.write(f"mkdir -p {rep_dir}\n\n")

            # Merge replicate BED files
            tmp_dir    = f"{work_dir}/02.results/01.preprocessing"
            id_paths   = [f"{tmp_dir}/{sid}/{sid}.sorted.bed" for sid in sample_ids]
            merged_bed = f"{rep_dir}/{rep_name}.rep.merged.bed"
            W2.write(f"cat {' '.join(id_paths)} > {merged_bed}\n" if len(id_paths) > 1
                     else f"ln -s {id_paths[0]} {merged_bed}\n")

            W2.write(f"\nml pybedtools/0.10.0-foss-2023b\n")

            # ── 77 bp window ──────────────────────────────────────────────────
            slop_bed    = f"{rep_dir}/{rep_name}.slop.bed"
            range77_bed = f"{rep_dir}/{rep_name}.range77.sorted.bed"
            bedgraph77  = f"{rep_dir}/{rep_name}.combined77_sorted.bedgraph"
            W2.write(f"bedtools slop -i {merged_bed} -g {chrom_sizes_slop} -b 2 -s > {slop_bed}\n")
            W2.write(f"python3 {script_path}fix_bed.py --input {slop_bed} --side l --length 77 --strandFlag True --output {range77_bed}\n")
            W2.write(f"bedtools genomecov -i {range77_bed} -bg -g {chrom_sizes_slop} > {bedgraph77}\n")

            # ── 10 bp window (FASTA + abundance) ─────────────────────────────
            slop6_bed   = f"{rep_dir}/{rep_name}.slop6.bed"
            range10_bed = f"{rep_dir}/{rep_name}.range10.sorted.bed"
            fa          = f"{rep_dir}/{rep_name}.fa"
            up_fa       = f"{rep_dir}/{rep_name}_upper.fa"
            nt_abun     = f"{rep_dir}/{rep_name}_nt_abun.txt"
            kmer_abun   = f"{rep_dir}/{rep_name}_kmer_abun.txt"
            W2.write(f"\nbedtools slop -i {merged_bed} -g {chrom_sizes_slop} -b 6 -s > {slop6_bed}\n")
            W2.write(f"python3 {script_path}fix_bed.py --input {slop6_bed} --side l --length 10 --strandFlag True --output {range10_bed}\n")
            W2.write(f"bedtools getfasta -fi {ref} -bed {range10_bed} -fo {fa} -s\n")
            W2.write(f"cat {fa} | tr 'a-z' 'A-Z' > {up_fa}\n")
            W2.write(f"python3 {script_path}fa2nt_abundance.py --input {up_fa} --output {nt_abun} --order TCGA --percentage\n")
            W2.write(f"python3 {script_path}kmer_abundance.py  --input {up_fa} --output {kmer_abun} --kmer 2 --percentage\n")

            # ── 2 bp window (FASTA + abundance) ──────────────────────────────
            slop2_bed   = f"{rep_dir}/{rep_name}.slop2.bed"
            range2_bed  = f"{rep_dir}/{rep_name}.range2.sorted.bed"
            bedgraph2   = f"{rep_dir}/{rep_name}.combined2_sorted.bedgraph"
            fa2         = f"{rep_dir}/{rep_name}.2fa"
            up2_fa      = f"{rep_dir}/{rep_name}_upper2.fa"
            nt_abun2    = f"{rep_dir}/{rep_name}_nt_abun2.txt"
            kmer_abun2  = f"{rep_dir}/{rep_name}_kmer_abun2.txt"
            W2.write(f"\nbedtools slop -i {merged_bed} -g {chrom_sizes_slop} -b 2 -s > {slop2_bed}\n")
            W2.write(f"python3 {script_path}fix_bed.py --input {slop2_bed} --side l --length 2 --strandFlag True --output {range2_bed}\n")
            W2.write(f"bedtools genomecov -i {range2_bed} -bg -g {chrom_sizes_slop} > {bedgraph2}\n")
            W2.write(f"bedtools getfasta -fi {ref} -bed {range2_bed} -fo {fa2} -s\n")
            W2.write(f"cat {fa2} | tr 'a-z' 'A-Z' > {up2_fa}\n")
            W2.write(f"python3 {script_path}fa2nt_abundance.py --input {up2_fa} --output {nt_abun2} --order TCGA --percentage\n")
            W2.write(f"python3 {script_path}kmer_abundance.py  --input {up2_fa} --output {kmer_abun2} --kmer 2 --percentage\n")

# SLURM array submission script for Step 2
with open(f"{work_dir}/01.scripts/02.replicates_merging/sub_array.sh", "w") as O2:
    O2.write('\n'.join([
        '#!/bin/bash',
        '#SBATCH --job-name=Damageseq_02.replicates_merging_array',
        '#SBATCH --partition=batch',
        '#SBATCH --ntasks=1',
        '#SBATCH --time=120:00:00',
        '#SBATCH --mem=64G',
        '#SBATCH --output=%x_%j.out',
        '#SBATCH --mail-type=END,FAIL',
        f'#SBATCH --array=1-{len(merged_dict)}',
        'file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" script.list)',
        'sh $file',
    ]) + "\n")
