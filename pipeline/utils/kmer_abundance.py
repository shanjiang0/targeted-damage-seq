#!/bin/python

import itertools
import sys
import argparse

def nucleotideAbundanceDict2percentageDict(theDict):
    newDict = {}
    for key in theDict.keys():
        newDict[key] = {}
        positionTotalCount = 0
        for letter in theDict[key].keys():
            positionTotalCount += theDict[key][letter]
        for letter in theDict[key].keys():
            newDict[key][letter] = float(theDict[key][letter])/positionTotalCount
    return newDict

def getMaxSeqLength(input_file):

    sequence = ''
    maxLength = 0
    filein = open(input_file, 'r')
    for line in filein:
        if line.startswith('>'):
            maxLength = max(maxLength, len(sequence))
            sequence = ''
        else:
            sequence += line.strip()
    else:
        maxLength = max(maxLength,len(sequence))
    filein.close()
    return maxLength

def getKmerAbundance(input_file,kmer, firstNletters=None):
    nucleotides = 'ATGC'
    kmer=int(kmer)
    subseqTupleList = itertools.product(nucleotides, repeat=int(kmer))
    subseqList = sorted(list(subseqTupleList))

    def subseqHasValidNucleotides(subseq):
        allowedLetters = 'ATCGatcg'
        for letter in subseq:
            if letter not in allowedLetters:
                return False
        return True

    def addSeqToPosDict(sequence, theDict):
        for i in range(min(firstNletters,len(sequence)) - (kmer - 1)):
            subseq = sequence[i : i + kmer]
            if subseqHasValidNucleotides(subseq):
                     theDict[i + 1][subseq.upper()] += 1

    if firstNletters == None:
        firstNletters = getMaxSeqLength(input_file)
    else:
        firstNletters = int(firstNletters)

    positionDict = {}
    for i in range(1, firstNletters + 1 - (kmer - 1)):
        positionDict[i] = {}
        for subseq in subseqList:
            positionDict[i][''.join(subseq)] = 0


    filein = open(input_file, 'r')
    seq = header = ''
    for line in filein:
        if line.startswith('>'):
            header = line.strip()
            addSeqToPosDict(seq, positionDict)
            seq = ''
        else:
            seq += line.strip()
    else:
        addSeqToPosDict(seq,positionDict)

    return positionDict
                

def writeKmerAbundanceTable(kmerAbundanceDict, output, percentage=False, append=False):
    if append:
        writeMode = 'a'
    else:
        writeMode = 'w'
    if percentage:
        kmerAbundanceDict = nucleotideAbundanceDict2percentageDict(kmerAbundanceDict)
    separator = '\t'
    out = open(output, writeMode)
    positions = kmerAbundanceDict.keys()
    for position in positions:
        out.write(separator + str(position))
    out.write("\n")
    for kmer in sorted(kmerAbundanceDict[1].keys()):
        out.write(kmer)
        for position in positions:
            out.write(separator + str(kmerAbundanceDict[position][kmer]))
        out.write("\n")
    return 1

def main():
    parser = argparse.ArgumentParser(description='abundance of kmer')
    parser.add_argument('--input', required = True, help='input_file')
    parser.add_argument('--kmer', required = True, help='')
    parser.add_argument('--percentage', action='store_true', help='Write percentages instead if actual counts')
    parser.add_argument('--output', required = True, help='output_file')

    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    kmer = args.kmer
    percentage = args.percentage
    positionDict = getKmerAbundance(input_file, kmer)
    writeKmerAbundanceTable(positionDict, output_file, percentage)

if __name__ == '__main__' :
    main()
