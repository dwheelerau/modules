#!/usr/bin/python

def gene_ontology(ontology_db):
    '''
    parses go term db and returns dic [GO:######]="info\tinfo\tinfo\t"
    found at copy at:
    http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/ontology/gene_ontology_edit.obo
    ...
    [Term]
    id: GO:0000001
    name: mitochondrion inheritance
    namespace: biological_process
    def: "The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]
    synonym: "mitochondrial inheritance" EXACT []
    is_a: GO:0048308 ! organelle inheritance
    is_a: GO:0048311 ! mitochondrion distribution

    [Term]
    id: GO:0000002
    name: mitochondrial genome maintenance
    namespace: biological_process
    def: "The maintenance of the structure and integrity of the mitochondrial genome; includes replication and segregation of the mitochondrial chromosome." [GOC:ai, GOC:vw]
    is_a: GO:0007005 ! mitochondrion organization
    '''
    #sadistic bastard!!! But works and fast!
    ontol_db = open(ontology_db).read()
    go_dic = {}
    data = ontol_db.split("[Term]")
    for go in data:
        info = go.strip().split("\n")
        #search for alt id
        alt_id = []
        for part in info:
            if part.find("alt_id:")!=-1:
                alt_id.append(part[part.find("GO:"):])    
        go_term = info.pop(0)[4:]
        go_info = "\t".join(info)
        go_dic[go_term]=go_info
        if alt_id:
            for alt in alt_id:
                go_dic[alt]=go_info
    return go_dic

        
        
