#!/usr/bin/python
'''
a group of tools that basically wraps up functions from the Biopython gff
module which is currently being developed and not part of the current biopython
release. See http://biopython.org/wiki/GFF_Parsing
It supports GFF 1,2,3!
'''
from BCBio.GFF import GFFParser
from sequence_utils import make_rc_string
from Bio import SeqIO


def get_feature_cord(gff_file,user_feature="gene"):
    """returns st, stop of a feature with orintation as tuple
    in a list to account for CDS
     in a dictionary with gene id as key:[(st,stop,orin)]
     freature can be one of 'protein','gene','mRNA','CDS','exon'
     user_feature="gene"
     a_dict['FBgn0031208']= [(7528, 9484, 1, '2L')]
     user_feature = "mRNA".
     a_dict['FBgn0031208']=[(7528, 9484, 1, 'FBtr0300689', '2L'),
      (7528, 9484, 1, 'FBtr0300690', '2L'),
       (7528, 9484, 1, 'FBtr0330654', '2L')]
    user_feature = "CDS"
    a_dict['FBgn0031208']= [([(7679, 8116), (8192, 8610)], 1, 'FBtr0300689', '2L'),
      ([(7679, 8116), (8192, 8589), (8667, 9276)], 1, 'FBtr0300690', '2L'),
        ([(7679, 8116), (8228, 8610)], 1, 'FBtr0330654', '2L')]
     """
    limit_info = dict(gff_type = ['protein','gene','mRNA','CDS','exon'])
    feature_dict = {}
    parser = GFFParser()
    in_handle = open(gff_file)
    for rec in parser.parse(in_handle,limit_info=limit_info):
        rec_id = rec.id
        for feat in rec.features:
            if feat.type == "gene":
                gene_id = feat.id
                if user_feature == "gene":
                    #gene_id = feat.id
                    assert gene_id not in feature_dict
                    feature_dict[gene_id]=[(feat.location.start.position,
                        feat.location.end.position,feat.strand,rec_id)]
                else:
                    for sub in feat.sub_features:
                        if sub.type == "mRNA":
                            if user_feature == "mRNA":
                                info = (sub.location.start.position,
                                        sub.location.end.position,sub.strand,sub.id,rec_id)
                                if gene_id in feature_dict:
                                    feature_dict[gene_id].append(info)
                                else:
                                    feature_dict[gene_id] = [info]
                            else:
                                codons = []
                                for sub_sub in sub.sub_features:
                                    if sub_sub.type == "CDS":
                                        st = sub_sub.location.start.position
                                        end = sub_sub.location.end.position
                                        codons.append((st,end))
                                info = (codons,sub.strand,sub.id,rec_id)
                                if gene_id in feature_dict:
                                    feature_dict[gene_id].append(info)
                                else:
                                    feature_dict[gene_id] = [info]
    in_handle.close()
    return feature_dict

def cds_dic(infile, GENBANK_COUNTING=None):
    '''inifle GFF record v1,2 or 3. Genbank counting (=1) will add a 1 to the first cord
    so that it uses DNA counting, default will return cords that can
    then be directly used to splice. Returns a dictionary using mRNA id
    ['mRNAid":(starnd,[list of codons as tuple pairs(st,end),(st,ed)]
    ,coding_seq (rev if req)]
    targets limits search to certain features, thus saving momory
     if sub.id == 
    biopython = [(271849, 272110)]
    genebank =  271850	- 272110
    >>>rec.seq.:tostring(biopythonstcord:biopython end cord]
    ATGGCAGCTGAGCAGTCAGTGTTGCTCGTCCAGCCTCAGACACAGCTGGTGAAGGTGGAGTCTCATCGCAGTTTTGTTTG
    CAAGCCAAAGAAGAAGCGCCTCCGCAAGGTGCTCAGCGAGAAGGAGAAATATTACAGACACCGTCGTGACTTCGAGCAGCGCATGGAGAAGCGTCTGGCCGGGATCGGCACAGTGCTGGCTAGGGTGATAGGTCATATATTAGTGCTGCAGTGTGTGTGGTCATGTGAAGCCGCTGTTTGA        
    for neg stand
    #Biopython = FBtr0118609 
    [(2794, 3160)]
    genbank
     2795	3160
    to collect would join exons then reverse compliment'''
    limit_info = dict(gff_type = ["protein", "gene", "mRNA", "CDS", "exon"],)
        #gff_type = ["protein", "gene", "mRNA", "CDS", "exon"],)
        #gff_id = ["scaffold_12944"],
        #)
    codon_dic = {}
    parser = GFFParser()
    in_handle = open(infile)
    for rec in parser.parse(in_handle, limit_info=limit_info):
        rec_seq = rec.seq.tostring()
        rec_id = rec.id
        for feature in rec.features:
            if feature.type == "gene":
                for sub in feature.sub_features:
                    if sub.type == "mRNA":
                        name = sub.id
                        codon_list = []
                        strand = sub.strand
                        seq_bits = []
                        for sub_sub in sub.sub_features:
                            if sub_sub.type == "CDS":
                                st = sub_sub.location.start.position
                                end = sub_sub.location.end.position
                                if GENBANK_COUNTING == None:
                                    codon_list.append((st,end))
                                    #WILL BREAK WITH GENBANK
                                    seq_bits.append(rec_seq[st:end])
                                    
                                else:
                                    #add 1 to st cord
                                    st = st + 1
                                    codon_list.append((st,end)) 
                        assert name not in codon_dic
                        #this currently only works for non-genebank type d
                        sequence = "".join(seq_bits)
                        if strand == 1:#type int
                            codon_dic[name] = [rec_id,strand,codon_list,sequence]
                        else:
                            assert strand == -1
                            #reverse compliment
                            codon_dic[name] = [rec_id,strand,codon_list,
                                            make_rc_string(sequence)]
                    
                            
    in_handle.close()
    return codon_dic

def scaffold_to_mRNA(gff_file):
    '''parses through and collects all the mRNA seq names on a scaffold
    returns a dic with scaffold_name:[list of mRNAs]'''
    limit_info = dict(
        gff_type = ["protein", "gene", "mRNA", "CDS", "exon"],)
        #gff_id = ["scaffold_12944"],
        #)
    mRNA_dic = {}
    parser = GFFParser()
    in_handle = open(gff_file)
    for rec in parser.parse(in_handle, limit_info=limit_info):
        #rec_seq = rec.seq.tostring()
        rec_id = rec.id
        mRNA_list = []
        for feature in rec.features:
            if feature.type == "gene":
                for sub in feature.sub_features:
                    if sub.type == "mRNA":
                        name = sub.id
                        if name not in mRNA_list:
                            mRNA_list.append(name)
        mRNA_dic[rec_id]=mRNA_list
    in_handle.close()
    return mRNA_dic

def extract_seq(gff_file,outfile):
    '''for gff with seq attached goes through and parses out to seq rec as
    fasta to a new file'''
    in_handle = open(gff_file)
    fasta_file = open(outfile,"w")
    parser = GFFParser()
    for rec in parser.parse(in_handle):#, limit_info=limit_info):
        #rec_seq = rec.seq.tostring()
        SeqIO.write(rec,fasta_file,"fasta")
    in_handle.close()
    fasta_file.close()
        
def get_feature_dict(infile,target='gene',
        limits=["protein", "gene", "mRNA", "CDS", "exon"]):


    """return dictionary of gene info from the gff file 'infile'
    dict= {gene_id:SeqIO.feature_record}
    target changes from gene to other features: ie
    mRNA, CDS, but feature tree will be more specific
    limits only touches those parts of the gff, may not always work
    """
    feature_dict = {}
    infile = open(infile)
    #limit_info = dict(gff_type = ["protein", "gene", "mRNA", "CDS", "exon"],)
    parser = GFFParser()
    for rec in parser.parse(infile,limit_info=dict(gff_type=limits,)):
        for feature in rec.features:
            if feature.type == target:
                feature_dict[feature.id]=feature
            else:
                pass
    return feature_dict

def extract_coding_seq(gff_file, fasta_file=None):
    #need to add version catch as version <3 dont have fasta file sometimes
    if fasta==None:
        codon_dic = cds_dic(gff_file)
        for mrna in codon_dic:
            pass
            
    else:
        pass#to be implemented
                    
                
    
