#!/usr/bin/python
def parse_gff(gff_file):
    ''' extracts cds and strand info from a gff file, returns this as a
    list and tuple['-', scffoldXX,('105230', '105408'), ('105719', '105840'),
    ('105918', '106033'), ('106223', '106270')] with key NV_name
    USAGE:strand=CDS_dic[name][0],coord=CDS_dic[name][3:],model =cds_dic[2]
    scaffold=CDS[1]***FIX:for multi transcripts just takes last one in gff file***'''
    gff = open(gff_file,"r")
    CDS_dic = {}
    for line in gff:
        line_info = line.split("\t")
        #start at id collect header info under gene name = key
        if line_info[8].find("ID=")!=-1:
            scaffold=line_info[0]
            model = line_info[1][7:]#removes version
            strand = line_info[6]
            gene_name = line_info[8][3:10]#for NV name are they all NVs?
            CDS_dic[gene_name]=[strand,scaffold,model]
        #not id line, then
        else:
            if line_info[2].find("CDS")!=-1:
                #update gene name just to be sure
                gene_name=line_info[8][7:14]
                #add exon information
                CDS_dic[gene_name].append((line_info[3],line_info[4]))#st/st
    gff.close()
    return CDS_dic

def slice_sequence(seq_string,CDS_tuple, strand="+"):
    """just slices a seq string based on cds info [and strand]"""
    mRNA="" #this should be converted to a seq object??
    counter = 0
    if strand!="+":
        for coord in CDS_tuple:
            counter+=1 
            st = coord[0]
            end = coord[1]
            #need to adjust for genetics counting versus pythonic counting
            if counter==1:#first segment which is really stop segment
                mRNA= mRNA+ seq_string[(int(st)-4):(int(end))]#Rev - add stop!
            else:
                mRNA= mRNA+ seq_string[(int(st)-1):(int(end))]
    else:
        counter = 0
        num_of_splice = len(CDS_tuple)
        for coord in CDS_tuple:
            counter +=1
            st = coord[0]
            end = coord[1]
            #need to adjust for genetics counting versus pythonic counting
            if counter!=num_of_splice:
                mRNA= mRNA+ seq_string[(int(st)-1):(int(end))]
            else:#last cds there need to catch stop codon
                mRNA= mRNA+ seq_string[(int(st)-1):(int(end)+3)]
    return mRNA
