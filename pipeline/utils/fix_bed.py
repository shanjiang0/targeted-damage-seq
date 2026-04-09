#!/usr/bin/env python
import sys
import argparse

def fixRange(input_file,side, length, output_file,strandFlag=True):
    I1 = open(input_file, 'r')
    W1 = open(output_file, 'w')
    for line in I1:
        content = bedLine2fixedRangedLine(line, side, length, strandFlag=True)
        W1.write(content + "\n")
    I1.close()
    W1.close()
    return

def bedLine2fixedRangedLine(bedLine, side, length,strandFlag=True):
    ll = bedLine.strip().split('\t')

    beg1 = int(ll[1])
    end1 = int(ll[2])
    #beg2 = int(ll[4])  # mate 2 start
    #end2 = int(ll[5])  # mate 2 end
    strand1 = ll[5]
    #strand2 = ll[9]
    #newBeg = 0
    #newEnd = 0
    #newBeg, newEnd = beg, end
    if (strandFlag == True and strand1 == '+') or (strandFlag == False):
        if side == 'l':
            newEnd1 = min(beg1 + length, end1)
            newBeg1 = beg1
        elif side == 'r':
            newBeg1 = max(end1 - length, beg1)
            newEnd1 = end1
    elif strandFlag == True and strand1 == '-':
        if side == 'l':
            newBeg1 = max(end1 - length, beg1)
            newEnd1 = end1
        elif side == 'r':
            newEnd1 = min(beg1 + length, end1)
            newBeg1 = beg1

    # For the second mate (mate 2), repeat the process
    #if (strandFlag == True and strand2 == '+') or (strandFlag == False):
    #    if side == 'l':
    #        newEnd2 = min(beg2 + length, end2)
    #        newBeg2 = beg2
    #    elif side == 'r':
    #        newBeg2 = max(end2 - length, beg2)
    #        newEnd2 = end2
    #elif strandFlag == True and strand2 == '-':
    #    if side == 'l':
    #        newBeg2 = max(end2 - length, beg2)
    #        newEnd2 = end2
    #    elif side == 'r':
    #        newEnd2 = min(beg2 + length, end2)
    #        newBeg2 = beg2

    if side != 'l' and side != 'r':
        sys.exit('Error: side must be either l or r. Exiting...')
    #print(f"Original: beg={beg}, end={end}, New: newBeg={newBeg}, newEnd={newEnd}")  # Debugging line
    ll[1] = str(newBeg1)
    ll[2] = str(newEnd1)
    #ll[4] = str(newBeg2)  # mate 2 new start
    #ll[5] = str(newEnd2)  # mate 2 new end
    return '\t'.join(ll)

def main():
    parser = argparse.ArgumentParser(description='fix fragment length if possible')
    parser.add_argument('--input', required = True, help='input_file')
    parser.add_argument('--side', required = True, help='side')
    parser.add_argument('--length', required =  True, help='sequence length of interest')
    parser.add_argument('--strandFlag', required = True, help='strand')
    parser.add_argument('--output', required = True, help='output_file')
    
    args = parser.parse_args()
    
    bedFile = args.input
    output_file = args.output
    side = args.side
    length = int(args.length)
#    strandFlag = args.strandFlag
    strandFlag = args.strandFlag.lower() == 'true'

    fixRange(bedFile,side, length, output_file,strandFlag)

if __name__ == "__main__" :
    main()

