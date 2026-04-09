#!/usr/bin/env python

import argparse
import re

# Function to convert nucleotide counts to percentages
def nucleotideAbundanceDict2percentageDict(theDict):
    newDict = {}
    for key in theDict.keys():
        newDict[key] = {}
        positionTotalCount = sum(theDict[key].values())  # Total counts at this position
        for letter in theDict[key]:
            newDict[key][letter] = float(theDict[key][letter]) / positionTotalCount
    return newDict

# Function to get the maximum sequence length in the FASTA file
def getMaxSeqLength(input_file):
    sequence = ''
    maxLength = 0
    with open(input_file, 'r') as filein:
        for line in filein:
            if line.startswith('>'):
                maxLength = max(maxLength, len(sequence))  # Update maxLength if new sequence is longer
                sequence = ''  # Reset sequence after encountering a new header line
            else:
                sequence += line.strip()
        maxLength = max(maxLength, len(sequence))  # Ensure the last sequence length is considered
    return maxLength

# Function to calculate nucleotide abundance (counts) at each position
def getNucleotideAbundance(input_file, firstNletters=None):

    if firstNletters is None:
        firstNletters = getMaxSeqLength(input_file)
    else:
        firstNletters = int(firstNletters)

    def letterIsValidNucleotide(letter):
        return letter.upper() in 'ATGC'  # Simplified nucleotide validation

    def addSeqToPosDict(sequence, theDict):
        for i in range(min(firstNletters, len(sequence))):
            letter = sequence[i]
            if letterIsValidNucleotide(letter):
                theDict[i + 1][letter.upper()] += 1

    positionDict = {i: {'A': 0, 'T': 0, 'G': 0, 'C': 0} for i in range(1, firstNletters + 1)}

    with open(input_file, 'r') as filein:
        seq = ''
        for line in filein:
            if line.startswith('>'):
                addSeqToPosDict(seq, positionDict)
                seq = ''  # Reset for the next sequence
            else:
                seq += line.strip()
        addSeqToPosDict(seq, positionDict)  # Process the last sequence

    return positionDict

# Function to write the nucleotide abundance (counts or percentages) to an output file
def writeNucleotideAbundanceTable(nucleotideAbundanceDict, output, nucleotideOrder='ATGC', percentage=False):
    if percentage:
        nucleotideAbundanceDict = nucleotideAbundanceDict2percentageDict(nucleotideAbundanceDict)
    
    nucleotides = nucleotideOrder
    with open(output, 'w') as out:
        out.write('\t' + '\t'.join(str(pos) for pos in nucleotideAbundanceDict.keys()) + '\n')
        for nucleotide in nucleotides:
            out.write(nucleotide)
            for pos in nucleotideAbundanceDict:
                out.write(f'\t{nucleotideAbundanceDict[pos][nucleotide]}')
            out.write('\n')

# Main function to handle argument parsing and calling appropriate functions
def main():
    parser = argparse.ArgumentParser(description='Calculate nucleotide abundance at each position in the sequence.')
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--order', required=True, help='Order of nucleotides, e.g., ATGC')
    parser.add_argument('--percentage', action='store_true', help='Output percentages instead of raw counts')
    parser.add_argument('--output', required=True, help='Output file')
    parser.add_argument('--firstNletters', type=int, help='Limit to the first N positions (optional)')
    
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    nucleotideOrder = args.order
    percentage = args.percentage
    firstNletters = args.firstNletters

    # Get the nucleotide abundance
    positionDict = getNucleotideAbundance(input_file, firstNletters)
    
    # Write the results to an output file
    writeNucleotideAbundanceTable(positionDict, output_file, nucleotideOrder, percentage)

if __name__ == '__main__':
    main()

