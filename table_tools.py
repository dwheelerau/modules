def bite(row,split_by="\t"):
    '''Strips new line only and splits (tab)a line args=any char'''
    return row.strip("\n").split(split_by)

def split_table(row,delineated_by):
    bits = row.split(delineated_by)
    return bits    

def strip_row(row,list_of_tags):
    for find in list_of_tags:
        row = row.strip(find)
    return row

def make_tab_dic(filename,delineated_by,row_key,remove_key="yes",list_of_rows="ALL",tup="NO"):
    '''all will take all line excpet key column, returns dic and error_list-[]'''
    file_handle = open(filename,"r")
    tab_dic = {}
    flag = 0
    repeat_list = []
    for row in file_handle:
        bits = split_table(row,delineated_by)
        key=bits[row_key]
        rows_to_add = []
        if list_of_rows=="ALL":
            columns = len(bits)
            rows_to_add = range(columns)
            if remove_key=="yes":
                rows_to_add.remove(row_key)#remove the key column)
        else:
            rows_to_add = list_of_rows[:]
        vals = []
        for targets in rows_to_add:
            vals.append(bits[targets])
        if tup!="NO":
            vals=tuple(vals)
        if key not in tab_dic:
            tab_dic[key]=vals            
        else:
            tab_dic[key]=tab_dic[key]+vals
            repeat_list.append(key)
    repeats = len(repeat_list)
    if repeats>0:
        #print "found %d keys with multiple vals, these were added together"%(repeats)
        #print "check these:"
        #print repeat_list
        return tab_dic,repeat_list
    else:
        return tab_dic,[]
    file_handle.close()
    

def rem_header(fileobj):
    '''
    Takes a fileobject with a header string that needs to be removed for parsing
    returns a tuple of file object (minus trailing \nwithout header and
    header string
    '''
    tmp = open("tmp.txt","w")
    flag = None
    header = ""
    for line in fileobj:
        if flag:
            if line!="\n":
                tmp.write(line)
        else:
            header = header + line
            flag = 1
    tmp.close()
    fileobj.close()
    newfile = open("tmp.txt")
    return (newfile,header)
            


    
