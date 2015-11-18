from blast_tools import strip_it
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os,sys
from file_scans import *
import re

def extract_sequences(seq_list,seq_file,output_file):

    '''extracts sequeces from the seq_file according to a list of seq
    names in seq_list. Sequences are saved in Fasta format in output_file.
    The fasta names in the header may contain additional information
    that was stored in the sequence object (.id)'''

    #need to say where the files are located windows problem commentout
    #path=os.getcwd()
    #full_path = path+'\\' + seq_file
    #out_path = path+'\\' + output_file

    #use function to make dictionary key=seq_name val=seqIO object
    seq_dic = collect_DNA_seq(seq_file)

    #output file
    out_handle = open(output_file,'w')#open(out_path,'w')
    
    for name in seq_list:
        #seq_dic[name].record.id = name
        try:
            out_handle.write(seq_dic[name].format("fasta"))
        except KeyError:
            print 'couldnt find sequence: ', name,\
                  'check for > or other format issues'
    #cleanup files
    out_handle.close()

def count_sequences(filename):
    '''
    counts fasta formated sequences such that returns an int containing
    the number of ">" ie sequences found in the file.
    '''
    file_handle=open(filename,'r')
    seq_counter = 0
    for line in file_handle:
        if line[0]=='>':
            seq_counter +=1
    file_handle.close()
    return seq_counter

def list_seq_names(filename):
    '''
    returns a list of cleaned sequence headers from a file of fasta sequences
    '''
    file_handle = open(filename)
    seq_list = []
    for line in file_handle:
        if line[0]=='>':
            seq_name = line.split(' ')[0]
            seq_name = strip_it(seq_name)
            seq_list.append(seq_name)
    print 'returned list is %d long'%(len(seq_list))
    return seq_list

def compare_lists(filename1,filename2):
    '''
    This function takes two files, makes a list of each and then
    and compares them, it reports on the number
    and names of sequences missing from either list, returns a list of missing
    sequences that can be extracted latter.
    '''

    seq_list1 = list_seq_names(filename1)
    seq_list2 = list_seq_names(filename2)
    seq_missing1 = [] #seq in list 1 not seen in list 2
    seq_missing2 = [] #seq in list 2 not seen in seq 1
    for names in seq_list1:
        if names not in seq_list2:
            seq_missing1.append(names)
    for names in seq_list2:
        if names not in seq_list1:
            seq_missing2.append(names)
    most_missing = []
    seq1_len,seq2_len = len(seq_missing1),len(seq_missing2)
    if seq1_len > seq2_len:
        most_missing = seq_missing1
        print 'lengths:file1(%d), file2(%d)'%(seq1_len,seq2_len)
    if seq1_len < seq2_len:
        most_missing = seq_missing2
        print 'lengths:file1(%d), file2(%d)'%(seq1_len,seq2_len)
    if seq1_len==seq2_len:
        print 'both the same length:file1(%d), file2(%d)'%(seq1_len,seq2_len)    

    return most_missing
        
def collect_seq_len(filename):

    '''requires "from Bio import SeqIO"'''
    
    file_handle = open(filename,'r')
    len_dic = {}

    for record in SeqIO.parse(file_handle,'fasta'):
        name = record.id
        seq = record
        if len(name.split(' '))>1:
               name = name.split(' ')[0].strip('>') #just fastsa header
        len_dic[name]=len(seq)
    file_handle.close()
    return len_dic

def collect_DNA_seq(fasta_file):
    '''returns a sequence dictionary containing all the dna sequences
    if a filename (fasta seqs) passed to this function. The key will be the seq
    name and the value will be a seqIO object from biopython. For example
    sequence could be written out:note inlclues fasta header'''

    seq_dic = SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta"))
    return seq_dic
    
def collect_raw_DNA(filename):
    '''like collect_DNA_seq but this time is a dicationary and the raw seq
    as a single string with no fasta markup'''
    file_handle = open(filename,'r')
    seq_dic={}
    for record in SeqIO.parse(file_handle,'fasta'):
        name = record.id
        if len(name.split(' '))>1:
               name = name.split(' ')[0].strip('>') #just fastsa header
        
        seq_dic[name] = record.seq.tostring()
    file_handle.close()
    file_handle.close()
    return seq_dic
    
def remove_stop_codons(dna_seq):#use with caution!
    '''removes typical stop codons from the ends of dna sequences. Takes a
    dna sequence as a string and returns a modified dna seq as a string'''
    stop_codons=["TAA","TAG","TGA"] #standard genetic code
    if dna_seq[-3:] in stop_codons:
        dna_seq=dna_seq[:-3]
    return dna_seq

def split_dna(dna_seq):
    '''breaks a dna seq string into a list of triplets (codons) ["atg","taa"]'''
    codons = []
    triplet = ""
    for nucleotide in dna_seq:
        triplet =triplet + nucleotide
        if len(triplet)==3:
            codons.append(triplet)
            triplet = ""
    return codons

def remove_third_positions(seq_file, outfile_name):
    '''first puts sequences in good format for parsing (removes carrage return
    etc before removing thrid codon positions. Will check that all seq are
    divisable by three ie coding regions'''
    temp_file = "temp.fna"
    inputFile=open(seq_file, 'r')
    outputFile=open(temp_file, 'w')
    for line in inputFile:
        if line[0] == '>':
            outputFile.write('\n' + line)
        else:
            outputFile.write(re.sub('[\d\s]','',line.upper()))
    inputFile.close()
    outputFile.close()
    
    #now check that all divisable by 3 ie coding regions
    len_dic = collect_seq_len(temp_file)
    for keys in len_dic:
        if len_dic[keys]%3!=0:
            print "%s is not divisable by three! Can't be coding DNA?"%(keys)
            
    #now remove thrid positions
    inputFile=open(temp_file, 'r')
    outputFile=open(outfile_name, 'w')
    for line in inputFile:
        counter=0
        #print fasta header
        if line[0]=="\n":
            pass
        elif line[0] == '>': 
            outputFile.write("\n"+line)
        else:
            for n in line: # from fasta.py each line == a full sequence
                counter+=1 
                if counter%3: # if remainder after division by three write to output
                    outputFile.write(n.upper())
    inputFile.close()
    outputFile.close()

def translate_file(seq_file,outfile):
    #I would not use this, use make_prot_record() below
    '''translates a file of fasta DNA sequences and writes them to file
    returns list of warning if an internal stop is encounted'''
    internal_stop=[]
    out_handle = open(outfile,"w")
    seq_dic = collect_DNA_seq(seq_file)
    #seq_dic = SeqIO_object_dic(seq_file)
    for names in seq_dic:
        name=seq_dic[names].id
        pep=Seq(str(seq_dic[names].seq.translate()))
        if pep[:-1].find("*")>=0:
            internal_stop.append(name)
        formated = ">"+name+"\n%s\n"%pep
        out_handle.write(formated)
    out_handle.close()
    return internal_stop
    
def fasta_formater(seq_string):
    '''a simple formater that prints 60 char per line for fasta
    takes a seq string as input ie str(seq) or seq.tostring()'''
    seq = ""
    step=60#format so that prints 60 char per line
    for i in range(0,len(seq_string),step):
        seq=seq+seq_string[i:i+step]+"\n"
    
    return seq  

def make_rc_record(record):
    """Returns a new SeqRecord with the reverse complement sequence.
    records = [make_rc_record(rec) for rec in SeqIO.parse("ls_orchid.fasta",
    "fasta")]    Now list comprehensions have a nice trick up their sleeves,
    you can add a conditional statement:    records = [make_rc_record(rec)
    for rec in SeqIO.parse("ls_orchid.fasta", "fasta") if len(rec)<700]
    """

    return SeqRecord(seq = record.seq.reverse_complement(),\
                     id = "rc_" + record.id, description = "reverse complement")

def make_rc_string(seq_string):
    """returns the rc of a sequence string as a string"""
    my_seq = Seq(seq_string)
    return my_seq.reverse_complement().tostring()

def fix_long_names(infile,outfile_seq,outfile_tab):
    """Blastdb and other programs hate long descriptive names, but these are
    handy sometimes, this scripts rips out all the crap and just calls them
    by there seqid, also outputs a translation table with all info displayed"""

    infile = open(infile)
    outfile = open(outfile_seq,"w")
    outfile_tab = open(outfile_tab,"w")

    for rec in SeqIO.parse(infile,"fasta"):
        name = rec.id
        dec = rec.description[len(name):]
        raw_seq = rec.seq.tostring()
        seq = fasta_formater(raw_seq)#converts to 60 char fasta
        new_fasta = ">%s\n%s\n"%(name,seq)
        outfile.write(new_fasta)
        table = "%s\t%s\n"%(name,dec)
        outfile_tab.write(table)
    infile.close()
    outfile.close()
    outfile_tab.close()

def make_record(seq_id,seq_str,description = ""):
    return SeqRecord(seq = seq_str, id = seq_id, description = description)

def make_index_dic(filename,seq_format="fasta"):
    user_dict = SeqIO.index(filename, seq_format)
    return user_dict

def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table).
    from Bio import SeqIO
    i.e
    proteins = (make_protein_record(nuc_rec) for nuc_rec in \
            SeqIO.parse("coding_sequences.fasta", "fasta"))
    SeqIO.write(proteins, "translations.fasta", "fasta")"""
    
    return SeqRecord(seq = nuc_record.seq.translate(cds=True), \
                     id = "trans_" + nuc_record.id, \
                     description = "translation of CDS, using default table")

def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table).
    i.e
    proteins = (make_protein_record(nuc_rec) for nuc_rec in \
            SeqIO.parse("coding_sequences.fasta", "fasta"))
    SeqIO.write(proteins, "translations.fasta", "fasta")"""
    return SeqRecord(seq = nuc_record.seq.translate(to_stop=True), \
                     id = "trans_" + nuc_record.id, \
                     description = "translation of CDS, using default table")
    
    
