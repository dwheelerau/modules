'''
a collection of functions for manipulating blast data
'''
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast.Applications import NcbitblastxCommandline

import sys,os
import subprocess
from linux_control import linux_controller

def strip_it(string):

    '''just strips out junk from a string, orig used strip method, now use replace'''
    
    clean = string
    clean = clean.replace(' ','')
    clean = clean.replace('>','')
    clean = clean.replace('\n','')
    clean = clean.replace('\t','')
    clean = clean.replace('+','')
    clean = clean.replace('-','')
    return clean

def find_seq_id(header):

    '''a function that takes header line and extracts
    just the name with no markups
    '''

    #for nomral type header wtih ' '
    header_bits = header.split(' ')
    header_name = header_bits[0]
    #cleanup the sequence removes '>'
    header_name = strip_it(header_name)
    return header_name
   
def parse_blast(blast_table_string):
    '''this will parse a blast table string and return all the details
    basically splits table by tabs and returns varibles for each.
    Returns a tuple (query(0),subject(1),%(2),overlap_len(3),mis_match(4)
    gaps (5),st_que_overlap(6), end q overlap(7), st_overlap_record(8),
    end_overlap(9),e-val(10),score(11)). 
    '''    
    line_bits = blast_table_string.split('\t')
        
    #collect all the info
    try:
        query_name,subject_name = line_bits[0],line_bits[1]
        percent = line_bits[2]
        overlap_len,miss_matches,gaps = line_bits[3],line_bits[4],line_bits[5]
        query_st,query_end = line_bits[6],line_bits[7]
        subject_st,subject_end = line_bits[8],line_bits[9]
        e_val,score = line_bits[10],strip_it(line_bits[11])
        tup_info = (query_name,subject_name,percent,overlap_len,miss_matches\
                ,gaps,query_st,query_end,subject_st,subject_end,e_val,score)
    except IndexError:
        tup_info = "posible problem with line %s"%(blast_table_string)
    return tup_info

def make_blastn_tab(query_file,db_file,outfile_name='blast_outfile.txt',e_cutoff=1e-5,max_return=1,max_aln=None,use_megablast=False,num_cores=1):
    '''Make a blastn table of DNA query and DNA database file with an e-val cutoff.
    Takes a quiry seq(s) as a file a blast db formated file using makeblastdb
    from the NCBI blast package, an e-value cutoff ie 1e-5 and an output file.
    max_return will be number returned. Output will be a blast table in the format[].
    The function essentially runs "blastn -query test.fna -db testdb
    -out test.xml -evalue 0.001 -outfmt 5" from the terminal!.
    ALT!!!! USE LINUX_CONTROL TO DO BLAST WITHOUT A DB IE TWO FAS FILES
    a = linux_controller("/home/dwheeler/software_tools/blast/bin/./blastn -query 'bait.txt' -subject 'test_db.txt' -out 'blast.out'") 
    a = None if works
    MASKING
    mask_cmd = blast_path+"./dustmasker -in '"+target_file+ "' -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out 'mask_info.asnb'"
    blastdb_cmd = blast_path+"./makeblastdb -in '"+target_file+"' -title 'target_db' -mask_data 'mask_info.asnb' -parse_seqids -out 'target_db' -dbtype nucl"
    
    '''

    #print warning if path different
    #print 'your blast path must be ~/software_tools/blast/bin/'
    cline = NcbiblastnCommandline(query=query_file, db=db_file, evalue=e_cutoff\
                                  , outfmt=6, out = outfile_name)

    #set up tasks
    if use_megablast:
        task = " -task megablast"
    else:
        task = " -task blastn"
    cpu = " -num_threads "+str(num_cores)    
    #number of alignments to return	
    if max_aln==None:
        other_opts = ' -max_target_seqs ' + str(max_return)
    else:
        #max_target_seqs is not compatable with num_desc or num_aln
        other_opts = " -num_alignments %s -num_descriptions %s"%(max_aln,max_aln)
    
    print 'processing this output...plase wait'
    command_sent = str(cline) + other_opts + task + cpu
    print command_sent
    
    return_code = subprocess.call(command_sent,shell=(sys.platform!="win32"))
    
    if return_code !=0:
        print 'sorry, we seem to have a problem!'
    else:
        print 'Done!....check %s when finished'%(outfile_name)
    #not sure how you know it is done


def make_blastx_tab(query_file,db_file,outfile_name='blast_outfile.txt',e_cutoff=1e-5, max_return=1,cpu=1):
    '''Make a blastn table of DNA query v's prot database file with an e-val cutoff.
    Takes a quiry seq(s) as a file a blast db formated file using makeblastdb
    from the NCBI blast package, an e-value cutoff ie 1e-5 and an output file.
    max_return will be number returned. Output will be a blast table in the format[].
    The function essentially runs "blastn -query test.fna -db testdb
    -out test.xml -evalue 0.001 -outfmt 5" from the terminal.
    '''
    #print warning if path different
    #print 'your blast path must be ~/software_tools/blast/bin/'   
    cpu = " -num_threads "+str(cpu)
    cline = NcbiblastxCommandline(query=query_file, db=db_file, evalue=e_cutoff\
                                  , outfmt=6, out = outfile_name)

    #set it up so it works
    
    max_return = str(max_return)
    other_opts = ('./',' -max_target_seqs'+' '+max_return)
    print 'processing this output...please wait'
    #command_sent = other_opts[0]+str(cline)+other_opts[1] + cpu
    command_sent = str(cline)+other_opts[1] + cpu
    print command_sent
    #file_path = '~/software_tools/blast/bin/'
    
    #return_code = subprocess.call(file_path+command_sent,shell=(sys.platform!="win32"))
    return_code = subprocess.call(command_sent,shell=(sys.platform!="win32"))
    
    if return_code !=0:
        print 'sorry, we seem to have a problem!'
    else:
        print 'Done!....check %s when finished'%(outfile_name)
    #not sure how you know it is done

def make_blastp_tab(query_file,db_file,outfile_name='blast_outfile.txt',e_cutoff=1e-5, max_return=1,cpu=1):
    '''Make a blastp table of prot v's prot database file with an e-val cutoff.
    Takes a quiry seq(s) as a file a blast db formated file using makeblastdb
    from the NCBI blast package, an e-value cutoff ie 1e-5 and an output file.
    max_return will be number returned. Output will be a blast table in the format[].
    The function essentially runs "blastn -query test.fna -db testdb
    -out test.xml -evalue 0.001 -outfmt 5" from the terminal.
    '''
    #print warning if path different
    #print 'your blast path must be ~/software_tools/blast/bin/'   
    cpu = " -num_threads "+str(cpu)

    cline = NcbiblastpCommandline(query=query_file, db=db_file, evalue=e_cutoff\
                                  , outfmt=6, out = outfile_name)

    #set it up so it works
    
    max_return = str(max_return)
    other_opts = ('./',' -max_target_seqs'+' '+max_return)
    print 'processing this output...please wait'
    #command_sent = other_opts[0]+str(cline)+other_opts[1] + cpu
    command_sent = str(cline)+other_opts[1] + cpu
    print command_sent
    #file_path = '~/software_tools/blast/bin/'
    
    #return_code = subprocess.call(file_path+command_sent,shell=(sys.platform!="win32"))
    return_code = subprocess.call(command_sent,shell=(sys.platform!="win32"))
    
    if return_code !=0:
        print 'sorry, we seem to have a problem!'
    else:
        print 'Done!....check %s when finished'%(outfile_name)
    #not sure how you know it is done

def make_tblastn_tab(query_file,db_file,outfile_name='blast_outfile.txt',
                     e_cutoff=1e-5, max_return=1,cpu=1):
    '''Make a blastp table of prot v's prot database file with an e-val cutoff.
    Takes a quiry seq(s) as a file a blast db formated file using makeblastdb
    from the NCBI blast package, an e-value cutoff ie 1e-5 and an output file.
    max_return will be number returned. Output will be a blast table in the format[].
    The function essentially runs "blastn -query test.fna -db testdb
    -out test.xml -evalue 0.001 -outfmt 5" from the terminal.
    '''
    #print warning if path different
    print 'your blast path must be ~/software_tools/blast/bin/'   
    cpu = " -num_threads "+str(cpu)

    cline = NcbitblastnCommandline(query=query_file, db=db_file, evalue=e_cutoff\
                                  , outfmt=6, out = outfile_name)

    #set it up so it works
    max_return = str(max_return)
    #other_opts = ('./',' -max_target_seqs'+' '+max_return)
    other_opts = (' -max_target_seqs'+' '+max_return)
    print 'processing this output...please wait'
    command_sent = other_opts[0]+str(cline)+other_opts[1] + cpu
    print command_sent
    #file_path = '~/software_tools/blast/bin/'
    
    return_code = subprocess.call(command_sent,shell=(sys.platform!="win32"))
    
    if return_code !=0:
        print 'sorry, we seem to have a problem!'
    else:
        print 'Done!....check %s when finished'%(outfile_name)
    #not sure how you know it is done

def blastn(query_file,subject_file,outfile_name='blast_outfile.txt',\
           e_cutoff="1e-5",outfmt="1",other_args=""):
    '''Blastn table two/more sequences use outfmt='6' for tab
    without having to make a database, suites small numbers of seqs 
    use outfmt=6 for table " -dust no" flag turrns off filter as other_args ie
    '''

    #print warning if path different
    #print 'your blast path must be ~/software_tools/blast/bin/'
    cline = " -query="+query_file+" -subject="+subject_file+" -outfmt="+outfmt +\
            " -out "+ outfile_name +" -evalue="+e_cutoff
    #set it up so it works
    #max_return = str(max_return)
    other_opts = ('./')#,' -max_target_seqs'+' '+max_return)
    #print 'processing this output...plase wait'
    command_sent = cline+other_args#+other_opts[1]
    file_path = '~/software_tools/blast/bin/./blastn'
    #print command_sent
    
    return_code = subprocess.call(file_path+command_sent,shell=(sys.platform!="win32"))
    
    if return_code !=0:
        return return_code
        #print 'sorry, we seem to have a problem!'
    else:
        pass #print 'Done!....check %s when finished'%(outfile_name)
    #not sure how you know it is done



def do_clustalw(fasta_seq_file):
    '''aligns multiple sequences in fasta format, saves file *.aln in
    your current working directory (containing the script'''
    cline = ClustalwCommandline("clustalw2", infile=fasta_seq_file)
    file_path = '~/software_tools/clustalw/./'+str(cline)
    print file_path
    return_code = subprocess.call(file_path,shell=(sys.platform!="win32"))
    if return_code !=0:
        #ERROR!
        return 1
    return None
    
def t_coffee_pairwise(fasta_seq_file):
    '''takes a infile of fasta sequences and writes a file to clustalw format'''
    cline = "t_coffee -infile "+fasta_seq_file
    file_path = "~/software_tools/T-COFFEE/bin/binaries/linux/./"+cline
    return_code = subprocess.call(file_path,shell=(sys.platform!="win32"))
    if return_code !=0:
        return 1
    return None

def find_multi_hits(blast_dic):
    '''this will take a blast dictionary and use anagram to find hits that
    appear twice, in other words it finds out what hits hit muliple contigs.
    The remainder are 1:1 matches. Returns a dictionary containing counts for
    the number of times each hit is pulled out of the database by the queries
    used. By accessing dic can get at 1:1 data and multiple data too
    '''
    anagram_dic = {} #dic to store all the data
    for names in blast_dic:        
        hit=blast_dic[names][0]
        if hit in anagram_dic:
            anagram_dic[hit]=anagram_dic[hit]+1
        else:
            anagram_dic[hit]=1
    return anagram_dic

#!/usr/bin/python

def recip_blast(sp1_blast_tab,sp2_blast_tab,outfile):
    '''takes two blast_tab_files with that only show top hit and compares
    them, returns a table file with only genes that are recip blast hits
    in both datasets, also prints number of excluded seqs from sp1_tab'''
    sp1_blast = open(sp1_blast_tab)
    sp2_blast = open(sp2_blast_tab)
    sp1_blast_hits, sp2_blast_hits = {}, {}

    #collect data in two dictionaries
    for line in sp1_blast:
        data = parse_blast(line)
        gene,hit = data[0],data[1]
        sp1_blast_hits[gene]=hit
    sp1_blast.close()    
    for line in sp2_blast:
        data = parse_blast(line)
        gene,hit = data[0],data[1]
        sp2_blast_hits[gene]=hit
    sp2_blast.close()

    #compare data in Nv hit to dros hit, are they reciprical?
    recip_hits = []
    dump=[]
    for gene in sp1_blast_hits: #shouldn't matter which dic
        sp1_hit = sp1_blast_hits[gene]
        try:
            sp2_hit = sp2_blast_hits[sp1_hit] #might not be any
        except KeyError:
            sp2_hit = None
        #gene = sp1## so compare gene to sp2 from dic
        if sp2_hit==gene:
            recip_hits.append(gene)
        else:
            #if ds not None then put hit in another list for checking
            if sp2_hit:
                dump.append(gene)
    #now go through Nv blast table and extract recip hits
    sp1_blast = open(sp1_blast_tab)
    recip_outfile = open(outfile,"w")
    for line in sp1_blast:
        gene= parse_blast(line)[0]
        if gene in recip_hits:
            recip_outfile.write(line)    
    sp1_blast.close()
    recip_outfile.close()
    #return the number of sequences  not recip in sp1 table
    return len(dump) 

def do_blat(dbase,qry,outfile="outfile",outform="psl",args=""):
    '''./blat 454Iso* cgb* -ooc=11.ooc output.psl
 psl - Default. See below other opt Tab separated format without actual sequence
-t=type     Database type.  Type is one of:
                 dna - DNA sequence
                 prot - protein sequence
                 dnax - DNA sequence translated in six frames to protein
               The default is dna
   -q=type     Query type.  Type is one of:
                 dna - DNA sequence
                 rna - RNA sequence
                 prot - protein sequence
                 dnax - DNA sequence translated in six frames to protein
                 rnax - DNA sequence translated in three frames to protein
               The default is dna
   -prot       Synonymous with -d=prot -q=prot
   -ooc=N.ooc  Use overused tile file N.ooc.  N should correspond to 
               the tileSize
   -tileSize=N sets the size of match that triggers an alignment.  
               Usually between 8 and 12
               Default is 11 for DNA and 5 for protein.
   -oneOff=N   If set to 1 this allows one mismatch in tile and still
               triggers an alignments.  Default is 0.
   -minMatch=N sets the number of tile matches.  Usually set from 2 to 4
               Default is 2 for nucleotide, 1 for protein.
   -minScore=N sets minimum score.  This is twice the matches minus the 
               mismatches minus some sort of gap penalty.  Default is 30
   -minIdentity=N Sets minimum sequence identity (in percent).  Default is
               90 for nucleotide searches, 25 for protein or translated
               protein searches.
   -maxGap=N   sets the size of maximum gap between tiles in a clump.  Usually
               set from 0 to 3.  Default is 2. Only relevent for minMatch > 1.
   -noHead     suppress .psl header (so it's just a tab-separated file)
   -makeOoc=N.ooc Make overused tile file
   -repMatch=N sets the number of repetitions of a tile allowed before
               it is marked as overused.  Typically this is 256 for tileSize
               12, 1024 for tile size 11, 4096 for tile size 10.
               Default is 1024.  Typically only comes into play with makeOoc
   -mask=type  Mask out repeats.  Alignments won't be started in masked region
               but may extend through it in nucleotide searches.  Masked areas
               are ignored entirely in protein or translated searches. Types are
                 lower - mask out lower cased sequence
                 upper - mask out upper cased sequence
                 out   - mask according to database.out RepeatMasker .out file
                 file.out - mask database according to RepeatMasker file.out
   -qMask=type Mask out repeats in query sequence.  Similar to -mask above but
               for query rather than target sequence.
   -minRepDivergence=NN - minimum percent divergence of repeats to allow 
               them to be unmasked.  Default is 15.  Only relevant for 
               masking using RepeatMasker .out files.
   -dots=N     Output dot every N sequences to show program's progress
   -trimT      Trim leading poly-T
   -noTrimA    Don't trim trailing poly-A
   -trimHardA  Remove poly-A tail from qSize as well as alignments in psl output
   -out=type   Controls output file format.  Type is one of:
                   psl - Default.  Tab separated format without actual sequence
                   pslx - Tab separated format with sequence
                   axt - blastz-associated axt format
                   maf - multiz-associated maf format
                   wublast - similar to wublast format
                   blast - similar to NCBI blast format
    '''
    command_str = "/home/dwheeler/software_tools/blat/./blat "\
                  +dbase+" "+qry+" -ooc=11.ooc " + outfile+"."+outform+" "+args
    print command_str
    os.system(command_str)

def make_tblastx_tab(query_file,db_file,outfile_name='blast_outfile.txt',
                     e_cutoff=1e-5, max_return=1):
    '''Make a blastp table of prot v's prot database file with an e-val cutoff.
    Takes a quiry seq(s) as a file a blast db formated file using makeblastdb
    from the NCBI blast package, an e-value cutoff ie 1e-5 and an output file.
    max_return will be number returned. Output will be a blast table in the format[].
    The function essentially runs "blastn -query test.fna -db testdb
    -out test.xml -evalue 0.001 -outfmt 5" from the terminal.
    '''
    #print warning if path different
    print 'your blast path must be ~/software_tools/blast/bin/'   

    cline = NcbitblastxCommandline(query=query_file, db=db_file, evalue=e_cutoff\
                                  , outfmt=6, out = outfile_name)

    #set it up so it works
    
    max_return = str(max_return)
    other_opts = ('./',' -max_target_seqs'+' '+max_return)
    print 'processing this output...please wait'
    command_sent = other_opts[0]+str(cline)+other_opts[1]
    print command_sent
    file_path = '~/software_tools/blast/bin/'
    
    return_code = subprocess.call(file_path+command_sent,shell=(sys.platform!="win32"))
    
    if return_code !=0:
        print 'sorry, we seem to have a problem!'
    else:
        print 'Done!....check %s when finished'%(outfile_name)
    #not sure how you know it is done

def make_nuc_db(seq,title = "target_db",folder=False,mask=True):
    #add code for making nuc and masking!
    blast_path = '~/software_tools/blast/bin/'
    #blastdb_cmd = blast_path+"./makeblastdb -in '"+seq+"' -title '"+title+"' -parse_seqids -out '"+title+"' -dbtype nucl"
    #linux_controller(blastdb_cmd)
    #subprocess.call([blastdb_cmd],shell=True)
    #if folder:
    #    linux_controller("mkdir "+folder+";mv "+title+".* "+folder)
    #add code for making nuc and masking!
    blastdb_cmd = blast_path+"./makeblastdb -in "+seq+" -title "+title+" -parse_seqids -out "+title+" -dbtype nucl"
    print blastdb_cmd
    #linux_controller(blastdb_cmd)
    subprocess.check_call(blastdb_cmd,shell=True)
    
    if folder:
        #linux_controller("mkdir "+folder+";mv "+title+".* "+folder)
        pass #not stable due to folder might allready exist
    
def make_pep_db(seq,title= "target_db",folder=False):
    #make a peptide db using makeblastdb save to folder if given
    blast_path = '~/software_tools/blast/bin/'
    blastdb_cmd = blast_path+"./makeblastdb -in "+seq+" -title "+title+" -parse_seqids -out "+title
    #linux_controller(blastdb_cmd)
    subprocess.check_call(blastdb_cmd,shell=True)
    if folder:
        #linux_controller("mkdir "+folder+";mv "+title+".* "+folder)
        pass #not stable
    
def extract_blastn_hits(seq_file,targets):
    #give seqfile, targets as blastn output will extract regions

    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file,"fasta"))

    outfile = open("hit_regions.fas","w")
    target_file = open(targets)
    
    for line in target_file:
        #print line
        bits = line.split("\t")
        target = bits[1]
        st_cord, end_cord = int(bits[8]),int(bits[9])
        #flip cord
        if st_cord>end_cord:
            tmp = st_cord
            st_cord = end_cord
            end_cord = tmp
        SeqIO.write(seq_dict[target][st_cord:end_cord],outfile,"fasta")

    outfile.close()
    target_file.close()    
