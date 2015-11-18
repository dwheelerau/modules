#!/usr/bin/env python

from Bio import SeqIO

#a group of functions for handling sequene slicing and cordinates
#ncbi counting starts at 1
#1234567
#ATGTTAA
#st = 1 end = 7 python would be 0:7

def ncbi_seq_slicer(ncbi_st,ncbi_end,seq_rec,rc=False):
    """return segment of sequence as seq record based on ncbi sequence cordinates
    takes int cordinates and biopython seqrecord. rc=True for reverse comp"""
    ncbi_st = int(ncbi_st)
    ncbi_end = int(ncbi_end)
    seq = seq_rec[ncbi_st-1:ncbi_end]
    
    if rc:
        seq_rc = seq.reverse_complement()
        seq_rc.id = seq.id
        seq_rc.description = seq.description
        return seq_rc
    else:
        return seq

if __name__ == "__main__":
    r = [rec for rec in SeqIO.parse('test.fasta','fasta')]
    a=ncbi_seq_slicer(1,5,r[0])




