#!/usr/bin/python
#####################################################
#   example.py - a program to ....                  #
#                                                   #
# Author: Dave Wheeler                              #
#                                                   #
# Purpose:                                          #
#                                                   #
# Usage:                                            #
#####################################################
#!/usr/bin/python
import sys
from Bio import SeqIO

def N50(numlist):
    """
    Abstract: Returns the N50 value of the passed list of numbers. 
    Usage:    N50(numlist)

    Based on the Broad Institute definition:
    https://www.broad.harvard.edu/crd/wiki/index.php/N50
    """
    numlist.sort()
    newlist = []
    for x in numlist :
        newlist += [x]*x
    
    # take the mean of the two middle elements if there are an even number
    # of elements.  otherwise, take the middle element
    #over 1/2 size of assembly in this contig size or bigger
    if len(newlist) % 2 == 0:
        medianpos = len(newlist)/2 
        return float(newlist[medianpos] + newlist[medianpos-1]) /2
    else:
        medianpos = len(newlist)/2
        return newlist[medianpos]

assert N50([2, 2, 2, 3, 3, 4, 8, 8]) == 6
assert N50([2, 2, 3, 4, 5, 6, 7]) == 5

if __name__=='__main__':
	lengths = []

	infile = open(sys.argv[1])

	for rec in SeqIO.parse(infile,'fasta'):
		lengths.append(len(rec))
	infile.close()

	print N50(lengths)
