#!/usr/bin/env python
import sys
import gzip

def ascii_to_phred(qual_str):
    """Convert ASCII-encoded quality scores to Phred scores."""
    return [ord(char) - 33 for char in qual_str]

def filter_fastq(read1_file, read2_file, output1_file, output2_file, quality_threshold=20):
    """Filter paired-end gzip-compressed FASTQ files based on Read 1 quality scores."""
    
    # Open input FASTQ files and output filtered files
    with gzip.open(read1_file, 'rt') as f1, gzip.open(read2_file, 'rt') as f2, \
         gzip.open(output1_file, 'wt') as out1, gzip.open(output2_file, 'wt') as out2:
        
        while True:
            # Read 4 lines (1 FASTQ record) from Read 1 and Read 2
            r1_header = f1.readline().strip()
            r1_seq = f1.readline().strip()
            r1_separator = f1.readline().strip()
            r1_quality = f1.readline().strip()

            r2_header = f2.readline().strip()
            r2_seq = f2.readline().strip()
            r2_separator = f2.readline().strip()
            r2_quality = f2.readline().strip()

            # Break if EOF is reached in any file
            if not r1_header or not r2_header:
                break

            # Convert quality string of Read 1 to Phred scores
            r1_phred_scores = ascii_to_phred(r1_quality)
            mean_phred_score_r1 = sum(r1_phred_scores) / len(r1_phred_scores)
            # Check if all quality scores in Read 1 are ≥ quality_threshold
            if mean_phred_score_r1 >= quality_threshold:
                # If Read 1 passes the quality check, write both Read 1 and Read 2 to output files
                out1.write(f"{r1_header}\n{r1_seq}\n{r1_separator}\n{r1_quality}\n")
                out2.write(f"{r2_header}\n{r2_seq}\n{r2_separator}\n{r2_quality}\n")

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 6:
        print("Usage: python script.py <read1_file> <read2_file> <output1_file> <output2_file> <quality_threshold>")
        sys.exit(1)

    # Get input/output files from command-line arguments
    read1_file = sys.argv[1]
    read2_file = sys.argv[2]
    output1_file = sys.argv[3]
    output2_file = sys.argv[4]
    quality_threshold = int(sys.argv[5])  # Quality threshold is passed as an integer

    # Run the filtering function
    filter_fastq(read1_file, read2_file, output1_file, output2_file, quality_threshold)

