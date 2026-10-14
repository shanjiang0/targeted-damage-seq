#!/bin/bash
set -euo pipefail

# Author      : Shan Jiang (sj31280@uga.edu)
# Affiliation : Li Lab, University of Georgia
# Repository  : https://github.com/shanjiang0/target-damage-seq
# Version     : 1.0.0
# License     : MIT (see LICENSE file)

# ── Project Setup Script for TARGET Damage-seq ────────────────────────────────
# Usage: bash build_project.sh
# Run this from your project working directory.
# Edit the configuration section below before running.

# Before you start this pipeline
# 1. Create the raw data folder:
#       mkdir -p 00.rawdata
# 2. Prepare your sample list and place it in 00.rawdata/.
#    The sample list is a strictly formatted file. See the example file in this repository.
# 3. Modify the configuration section below to match your local environment.
# 4. Run:
#       bash build_project.sh
# 5. All generated scripts will be placed in 01.scripts/
#    and all outputs will be written to 02.results/

# ── Configuration (edit these) ────────────────────────────────────────────────

METHOD="Targeted_Damageseq"
treatment="UV"
seq="PE"
workdir="$(pwd)"

# Reference paths — replace these with your local paths
bt_ref="/path/to/bowtie1_index/hg19"
ref_fa="/path/to/reference/hg19.fa"
chrom_sizes_slop="/path/to/reference/hg19.chrom.sizes"

# Pipeline script location — replace with your local path
pipeline_script="/path/to/target-damage-seq/pipeline/Targeted_Damageseq_PE_pipeline.py"

# ── Build directory structure ─────────────────────────────────────────────────

for i in "01.scripts" "02.results"; do
    for k in "01.preprocessing" "02.replicates_merging"; do
        mkdir -p "${workdir}/${i}/${k}"
    done
done
mkdir -p "${workdir}/00.rawdata"

echo "Directory structure created."

# ── Sample list ───────────────────────────────────────────────────────────────

sample_list="${workdir}/00.rawdata/targeted_${METHOD}_sample.list"

if [ ! -f "$sample_list" ]; then
    echo "ERROR: Sample list not found at $sample_list"
    echo "Please create it using the example format in targeted_Damageseq_sample.list.example"
    exit 1
fi

# ── Generate pipeline run script ──────────────────────────────────────────────

cd "${workdir}/01.scripts/01.preprocessing"

runsh="${workdir}/01.scripts/01.preprocessing/get_${METHOD}_${seq}_pipeline_run.sh"

cat > "$runsh" <<EOF
#!/bin/bash
python3 "$pipeline_script" "$sample_list" "$workdir" "$treatment" "$bt_ref" "$ref_fa" "$chrom_sizes_slop"
EOF

chmod +x "$runsh"
bash "$runsh"


echo "If generated successfully, submit jobs with:"
echo "  cd ${workdir}/01.scripts/01.preprocessing && sbatch sub_array.sh"
echo "  cd ${workdir}/01.scripts/02.replicates_merging && sbatch sub_array.sh"
