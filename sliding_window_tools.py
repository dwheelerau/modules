#!/usr/bin/python
#import filter_order_output

def each_cons_iter(listin,n): 
    """(overlapp'g)moving window,return n items at a time from listin thru iterator"""
    i=0
    while (i<len(listin)-n+1 ):
        yield listin[i:i+n]
        i+=1
        
def each_slice_iter(listin,n):    
    """non-overlapp'g slices n elements at a time from listin,returned thru iterator"""
    i=0;        len_listin=len(listin)
    while (i<len_listin ):
        yield listin[i:min(len_listin,i+n)]
        i+=n
        
def each_cons_lol(listin,n):    
    """moving window, return (list of lists) of n items at a time"""
    return [listin[i:i+n] for i in range(len(listin)-n+1)]

def each_slice_lol(listin,n):        
    """non-overlapp'g slices, return (list of lists) """
    len_listin=len(listin)
    return [listin[i:min(len_listin, i+n)] for i in range(0,len_listin,n)]

def get_pat(window,compare_dic):
    '''outputs a list of tuples, freq or '('1','0'),X) is 'X' of tuple
    [(('0', '1'), 1), (('1', '1'), 1), (('0', '1'), 1), (('1', '1'), 1), (('0', '1'), 1)]
    [(('0', '1'), 2), (('1', '1'), 1), (('0', '1'), 1), (('1', '1'), 1)]
    '''
    patterns = []
    for gene in window:
        patterns.append(compare_dic[gene])
    data = []
    count = 1
    for n in range(len(patterns)):
        try:
            if patterns[n]==patterns[n+1]:
                count+=1
            else:
                data.append((patterns[n],count))
                count=1
        except IndexError:#end of list, save last info
            data.append((patterns[n],count))
    return data

def sort_pattern(data,window):
    #sort data hightest to lowest by times found
    data_tmp = sorted(data,key=lambda pos: pos[1],reverse=True)
    common_pattern = data_tmp.pop(0)
    common_index = data.index(common_pattern) #index pos
    pos_count = 0
    counter=0
    for info in data:
        if counter==common_index:
            first_gene = pos_count
            last_gene = pos_count+info[1]
            genes = window[first_gene:last_gene]
            return((window,data,genes,common_pattern))
        else:
            counter+=1 #move to next index
            pos_count = pos_count+info[1]#add onto index tracker!

def checker(best,data_tmp,window):
    #is the next best the same score as the first best
    #FINISH! at the momment just prints, returns best accord to sort()
    next_val = data_tmp.pop(0)
    count = next_val[1]
    if count<best[1]:
        return None #all good is a smaller number
    else:#fix!!!
        if len(set(next_val[0]))>1:
            print "dead heat!:%s v's %s in %s"%(best,next_val,",".join(window))
            return None
        else:
            return None #all good is a same value number
        

def sort_pattern_opt(data,window,limit):
    #sort data hightest to lowest by times found
    data_tmp = sorted(data,key=lambda pos: pos[1],reverse=True)
    common_pattern = None

    #search through and make sure at least one difference (1,0) rather than(1,1)
    for pat in data: #loop for as many values
        best = data_tmp.pop(0)
        try:
            if len(set(best[0]))>1 and best[1]>=int(limit):
                common_pattern = best
                if len(data_tmp)>0:
                    check = checker(best,data_tmp,window)
                break #keep current val as set((1,1,1))=1 while set((1,2,2)=2
            else:
                pass
        except TypeError: #called if [none,none]
            pass
        
    if common_pattern:
        common_index = data.index(common_pattern) #index pos
        pos_count = 0
        counter=0
        for info in data:
            if counter==common_index:
                first_gene = pos_count
                last_gene = pos_count+info[1]
                genes = window[first_gene:last_gene]
                return((window,data,genes,common_pattern))
            else:
                counter+=1 #move to next index
                pos_count = pos_count+info[1]#add onto index tracker!
    else:
        return None
    
def make_array_dic(filename):
    #collect all the exrpesiion data for each gene
    array_dic = {}
    file_h = open(filename)
    for line in file_h:
        bits = line.strip().split("\t")
        assert bits[0] not in array_dic
        array_dic[bits[0]]=bits
    file_h.close()
    return array_dic

def sliding_win(ordered,get_list,out_file,table_file,gene_pos_dic,window_size,limit):
    #function that controls others
    array_dic = make_array_dic(table_file)
    tmp_file = open(out_file,"w")
    target_cols = [str(num) for num in get_list]#just for str formating
    tmp_file.write("1=expn 'ON',0=expn 'OFF'; comparisons (column from table(%s)=(%s)); window_size=%s;limit=%s genes\n"%\
                   (table_file,",".join(target_cols),window_size,limit))
    len_get_list = len(get_list)

    for n in range(len(ordered)):
        compare_dic = {}
        window = ordered[n:n+window_size]#[n:n+10]
        for name in window:
            if name not in array_dic:
                compare_dic[name]=[None]*len_get_list
            else:
                tmp = []
                for n in get_list:
                    tmp.append(array_dic[name][n])#access expn info from other dic
                compare_dic[name]=tuple(tmp)

        data = get_pat(window,compare_dic)
        results = sort_pattern_opt(data,window,limit)
        '''(['NV21560', 'NV10580', 'NV10581', 'NV10582', 'NV10583'], [(('1', '0'), 1)
        , (('1', '1'), 4)], ['NV10580', 'NV10581', 'NV10582', 'NV10583'],(('1', '1'), 4))
  
        '''
        #output genes,number cord genes found, pattern
        if results:
            positions = (gene_pos_dic[results[2][0]],gene_pos_dic[results[2][-1]])
            output="%s\t%s\t%s\t%s\n"%(results[3][0],results[3][1],positions\
                                       ,",".join(results[2]))
            tmp_file.write(output)
        #fail if window only has all on or all off, or fails limit test
        else:
            pass
    tmp_file.close()

#run this script
if __name__=="__main__":
    from config_file import *
    import gene_order #import ordered_gene_list
    gff_dic = CDS_functions.parse_gff(gff_file)
    ordered,cord_dic  = gene_order.ordered_gene_list(gff_dic,scaffold_target\
                                                     ="SCAFFOLD1")
    
    gene_pos_dic = {}
    for pos in cord_dic:#from gene_order.py
        gene_pos_dic[cord_dic[int(pos)]]=pos
    #for all scaffolds collect a different ordered list using ordered_gene_list       
    sliding_win(ordered,get_list,out_file,table_file,gene_pos_dic)
    filter_order_output.filter_order(out_file,"testing.txt")
    
    

