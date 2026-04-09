#!/usr/bin/env bash

# Input FASTQ files
read1_fastq=$1  # Input Read 1 FASTQ file (gzipped)
read2_fastq=$2  # Input Read 2 FASTQ file (gzipped)

# Output FASTQ files (gzipped)
output_read1=$3  # Output Read 1 FASTQ file
output_read2=$4  # Output Read 2 FASTQ file

# Use the process ID to create unique temporary file names
pid=$$
temp_read1="temp_read1_${pid}.txt"
temp_read2="temp_read2_${pid}.txt"
temp_combined="temp_combined_${pid}.txt"
deduplicated_combined="deduplicated_combined_${pid}.txt"

# Ensure output directories exist
mkdir -p "$(dirname "$output_read1")"
mkdir -p "$(dirname "$output_read2")"

# Step 1: Decompress both FASTQ files and combine Read 1 and Read 2 sequences into the same row
zcat "$read1_fastq" | paste - - - - > "$temp_read1"
zcat "$read2_fastq" | paste - - - - > "$temp_read2"
paste "$temp_read1" "$temp_read2" > "$temp_combined"

# Step 2: Sort and deduplicate based on Read 1 (col 2) and Read 2 (col 6) sequences
sort -u -k2,2 -k6,6 "$temp_combined" > "$deduplicated_combined"

# Step 3: Ensure output files are empty before writing (overwrite mode)
: > "$output_read1"
: > "$output_read2"

# Step 4: Split deduplicated combined file back into Read 1 and Read 2 FASTQ format
while read -r col1 col2 col3 col4 col5 col6 col7 col8; do
    # Write the Read 1 (columns 1-4) to the output file for Read 1
    echo -e "$col1\n$col2\n$col3\n$col4" >> "$output_read1"

    # Write the Read 2 (columns 5-8) to the output file for Read 2
    echo -e "$col5\n$col6\n$col7\n$col8" >> "$output_read2"
done < "$deduplicated_combined"

ml load tabix/0.2.6-GCCcore-13.3.0
# Step 5: Compress the output files using bgzip
bgzip -c "$output_read1" > "${output_read1}.gz"
bgzip -c "$output_read2" > "${output_read2}.gz"


# Cleanup temporary files
rm "$temp_read1" "$temp_read2" "$temp_combined" "$deduplicated_combined"

