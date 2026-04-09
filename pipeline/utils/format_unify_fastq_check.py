#!/usr/bin/env python
import sys
import gzip

def is_valid_sequence_line(sequence):
    """Check if the sequence line contains only A, T, C, G, N at the start."""
    return sequence and sequence[0] in 'ATCGN'

def is_valid_quality_line(quality):
    """Check if the quality line contains valid ASCII characters (ignoring full validation)."""
    return quality and all(33 <= ord(char) <= 126 for char in quality)

def clean_and_format_sequence(sequence):
    """Convert sequence to uppercase and replace non-ATCGN characters with N."""
    sequence = sequence.upper()  # Convert to uppercase
    sequence = ''.join([base if base in 'ATCGN' else 'N' for base in sequence])
    return sequence

def remove_tabs_and_spaces(line):
    """Remove all tabs and spaces from a line."""
    return line.replace('\t', '').replace(' ', '')

def process_fastq(read1_file, read2_file, output_read1_file, output_read2_file):
    """Process paired-end FASTQ files and output only valid blocks."""
    
    with gzip.open(read1_file, 'rt') as f1, gzip.open(read2_file, 'rt') as f2, \
         gzip.open(output_read1_file, 'wt') as out1, gzip.open(output_read2_file, 'wt') as out2:

        while True:
            # Process Read 1
            r1_header = f1.readline().strip()
            if not r1_header:  # EOF
                break

            # Process Read 2
            r2_header = f2.readline().strip()
            if not r2_header:  # EOF
                break

            # Step 1: Check that both Read 1 and Read 2 headers start with '@'
            if not (r1_header.startswith('@') and r2_header.startswith('@')):
                continue  # Skip to the next line until valid headers are found

            # Step 2: Process the sequence lines (only A, T, C, G, N allowed at the start)
            r1_seq = f1.readline().strip()
            r2_seq = f2.readline().strip()

            if not (is_valid_sequence_line(r1_seq) and is_valid_sequence_line(r2_seq)):
                continue  # Skip to the next block if sequences are invalid

            # Clean and format sequences
            r1_seq = clean_and_format_sequence(r1_seq)
            r2_seq = clean_and_format_sequence(r2_seq)

            # Step 3: Process the separator lines (must be '+' or '-')
            r1_separator = f1.readline().strip()
            r2_separator = f2.readline().strip()

            if not (r1_separator in ['+', '-'] and r2_separator in ['+', '-']):
                continue  # Skip to the next block if separators are invalid

            # Step 4: Process the quality score lines (ASCII characters at the start)
            r1_quality = f1.readline().strip()
            r2_quality = f2.readline().strip()

            if not (is_valid_quality_line(r1_quality) and is_valid_quality_line(r2_quality)):
                continue  # Skip to the next block if quality lines are invalid

            # Remove tabs and spaces from all lines
            r1_header = remove_tabs_and_spaces(r1_header)
            r1_seq = remove_tabs_and_spaces(r1_seq)
            r1_separator = remove_tabs_and_spaces(r1_separator)
            r1_quality = remove_tabs_and_spaces(r1_quality)

            r2_header = remove_tabs_and_spaces(r2_header)
            r2_seq = remove_tabs_and_spaces(r2_seq)
            r2_separator = remove_tabs_and_spaces(r2_separator)
            r2_quality = remove_tabs_and_spaces(r2_quality)

            # If the block is valid, write it to the output
            out1.write(f"{r1_header}\n{r1_seq}\n{r1_separator}\n{r1_quality}\n")
            out2.write(f"{r2_header}\n{r2_seq}\n{r2_separator}\n{r2_quality}\n")

def main():
    """Main entry point for the script."""
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 5:
        print("Usage: python script.py <read1_file.gz> <read2_file.gz> <output_read1_file.gz> <output_read2_file.gz>")
        sys.exit(1)

    # Get input/output files from command-line arguments
    read1_fastq = sys.argv[1]
    read2_fastq = sys.argv[2]
    output_read1_file = sys.argv[3]
    output_read2_file = sys.argv[4]

    # Process the FASTQ files and output only valid blocks
    process_fastq(read1_fastq, read2_fastq, output_read1_file, output_read2_file)

if __name__ == "__main__":
    main()

