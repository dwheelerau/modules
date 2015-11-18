#!/usr/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
from sequence_utils import *
import sys
import subprocess

def check_seq(pep_dic,dna_dic):
    '''compares a conceptual translation to the peptide sequence in the aligment
    to make sure that the dna sequence given actually encodes the protein under
    the same name'''
    translation_dic = {}
    #make a dic of conceptually translated sequences
    missing_list = []
    #crashing here because of gaps?
    for dna_names in dna_dic:        
        translation_dic[dna_names]=str(dna_dic[dna_names].seq.translate())
    for pep_names in pep_dic:
        errors=[]
        error_file ="error_log.txt"
        pep_seq = pep_dic[pep_names].seq.tostring().replace("-","")
        #check for name miss-matches and report then crash if found
        try:
            dna_translation = translation_dic[pep_names]
        #catch sequence protein seq that dont have a DNA partner
        except KeyError:
            errors.append(pep_names)
        if len(errors)>0:
            error_handle = open(error_file,"w")
            errors_message = "".join(errors)
            error_handle.write(errors_message)
            error_handle.close()
            print "DNA and protein seq names don't match, check %s"%(error_file)
            print "The program will now 'crash'!"
            exit() #crash the program
        #check that the two sequences are the same!
        if pep_seq.upper().strip("*") != dna_translation.strip("*"):
            print pep_seq,dna_translation
            missing_list.append(pep_names)
    if len(missing_list) == 0:
        return None #success!
    else:
        return missing_list #problems with AA sequence itself
    
def find_gap(AA_seq):
    '''make a list of gap positions'''
    counter = 0
    position_list = []
    for aa in AA_seq:
        if aa == "-":
            position_list.append(counter)
        counter+=1
    return position_list

def insert_gaps(aa_seq,dna_seq):
    '''split dna into codons then insert gaps based on gap list. By using
    codons 1 AA postion = 1 DNA (codon position) and you can map the sequences
    to each other. Note that no codon/AA matching is done, that is why I pass the
    data to the check_seq function first'''
    codons=split_dna(dna_seq)
    gap_list = find_gap(aa_seq)
    for numbers in gap_list:
        codons.insert(numbers,"---")
    edited_dna = "".join(codons)
    gapped_dna = remove_stop_codons(edited_dna)
    return gapped_dna


def gap_a_lator(prot_aln_dic,dna_dic):

    '''takes protein alignments as a dictionary and their corresponding
    open reading frame DNA (with the same name!) and aligns the DNA based
    on the protein sequence alignment. This will perform checks and report
    errors if anything funny is going on.
    '''
    
    missing_list = check_seq(prot_aln_dic,dna_dic)
    if missing_list != None:
        print "problems found in:",missing_list
    for name in prot_aln_dic:
        pep_seq = prot_aln_dic[name].seq.tostring().upper()
        dna_seq = dna_dic[name].seq.tostring().upper()
        gapped_dna = insert_gaps(pep_seq,dna_seq)
        dna_dic[name].seq = Seq(gapped_dna)
    outfile_name = "gapped.fna"
    outfile_handle = open(outfile_name,"w")
    for names in dna_dic:
        outfile_handle.write(dna_dic[names].format("fasta"))
    outfile_handle.close()

def pairwise_codeml():

    '''takes a fasta aligned file string and returns the dn/ds value as a string
    the control file needs to be in the same directory as the infile
    **requires that each file to be processed
    be called "infile" or what ever is in the the ctr file that must be with
    the sequences!! So open each file and save it as infile in a loop****'''
    
    paml_path = "/home/dwheeler/software_tools/paml44/bin/"
    codeml_command = paml_path+"./codeml.exe > checkfile.fas" #temp error file
    return_code = subprocess.call(codeml_command,shell=(sys.platform!="win32"))
    if return_code !=0:
        return None
    else:
        output_handle = open("rst").read()
        data = output_handle.find("Paras")+17 #catch data with dirty hack!
        return output_handle[data:].strip()
            
def pairwise_baseml():

    '''takes a fasta aligned file string and returns the divergencevalue as a
    string
    the control file needs to be in the same directory as the infile
    **requires that each file to be processed
    be called "infile" or what ever is in the the ctr file that must be with
    the sequences!! So open each file and save it as infile in a loop****'''
    
    paml_path = "/home/dwheeler/software_tools/paml44/bin/"
    baseml_command = paml_path+"./baseml.exe > checkfile.fas" #temp error file
    return_code = subprocess.call(baseml_command,shell=(sys.platform!="win32"))
    if return_code !=0:
        return None
    else:
        output_handle = open("rst").read()
        data = output_handle.find("Paras")+17 #catch data with dirty hack!
        return output_handle[data:].strip()
            
