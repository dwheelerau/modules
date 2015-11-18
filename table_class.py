#!/usr/bin/python
'''
a simple class for dealing with tables
'''

class Table:
    '''a basic table class of rows and columns'''
    def __init__(self,row):
        self.row = row.strip("\n")

    def check_line(self,separated_by):
        '''rules that could break program ie empty lines, returns none if probs'''
        if self.row.find(separated_by)==-1:
            return None
        if self.row =="\n":
            return None
        if self.row ==None:
            return None
        else:
            return 1

    def split_row(self,separated_by="\t"):
        return row.split(separated_by)
        
class Table_split(Table):
    '''a class that splits a table into two objects'''
    def split_for_dic(self,break_point,separate_by="\t"):
        '''returns strings with key and values separated as specif by user'''
        self.separate_by = separate_by
        key,val= self.split_list(break_point,self.separate_by)
        key,val = "".join(key),self.separate_by.join(val)
        return key,val
    
    def split_list(self,break_point, separate_by):
        '''specify column position to split'''
        bits = self.row.split(separate_by)
        start = bits[:break_point]
        end = bits[break_point:]
        return start, end
        
class Table_dic:
    '''returns dic with key as first value, use Yes to skip first line'''
    def __init__(self,filename,delineated_by,skip_first_line="No"):
        self.filename_handle = open(filename,"r")
        self.delineated_by = delineated_by
        self.tab_dic={}
        self.skip_first_line = skip_first_line

    def make_dic(self):
        first_line_test = 0
        if self.skip_first_line=="No":
            first_line_test=1
        for line in self.filename_handle:
            if first_line_test==0:
                pass
            else:
                row = Table_split(line)
                test = row.check_line(self.delineated_by)
                if test:
                    #hard wired change - fix make more general
                    key,val=row.split_for_dic(1,self.delineated_by)
                    self.tab_dic[key]=val
            first_line_test=1
        self.filename_handle.close()
        return self.tab_dic


class Specific_split(Table):
    '''specify a location (human format ie column 1 = 1) for a key and a key
    value par will be returned. The remaining row will be spliced back together
    '''
    def __init__(self,row,key_column, split_by):
        Table.__init__(self,row)
        self.key_column = key_column
        self.split_by = split_by

    def extract_key_val(self):
        try:
            self.key = self.row.split(self.split_by)[self.key_column-1]
            self.val = self.row.split(self.split_by)
            self.val.remove(self.key)
            return self.key, self.val
        except IndexError:
            return None,None

class Make_seq_dict:
    '''returns a dictionary with user option for key value, key is remvoed from
    "val" section'''
    def __init__(self, filename, key_column=1,split_by="\t"):
        self.filename = filename
        self.key_column = key_column
        self.split_by = split_by

    def make_new_dic(self):
        self.file_handle = open(self.filename)
        self.new_dic = {}
        for line in self.file_handle:
            data = Specific_split(line,self.key_column,self.split_by)
            key,val = data.extract_key_val()
            if key:
                self.new_dic[key]=val
        self.file_handle.close()
        return self.new_dic
            

class Join_tables:
    '''a class that joins tables, looks up columns and ids them automatically
    based on a list of common atrobutes: REQ ........

    def __init__(self,key_file,filename1,filename2,split_by="\t"):
        self.key_file=key_file
        self.filename1 = filename1
        self.filename2 = filename2
        self.handle1 = open(self.filename1)
        self.handle2 = open(self.filename2)
        self.split_by=split_by
        self.skip_line = 3
        self.count = 0
        #create list containng 100 colmns of data to check val tpyes
        self_list1 = create_short_row_list(self.handle1)
        self_list2 = create_short_row_list(self.handle2)
        list1_len = len(self_list1[1]) #skip first line just incase
        list2_len = len(self_list2[1])
        for n in range(list1_len):
            for val in list1[n]:
                pass
                
    
        #find row1

    def create_short_row_list(self,handle):
        self.counter=0
        self.handle = handle
        the_list = []
        for line in self.handle:
            the_list.append([bits for bits in line.split("\t")])
            self.counter+=1
            if self.counter>100:
                break
        return the_listdef find_good_keys(self,row):
        self.row = row
        self.key_handle = open(self.key_file)
        for line in self.key_handle:
            for bit in self.row:
                for bits in line:
                    if bit ==bits
            

    def find_good_key(self,row):#want to search through:
        #start with NV number
        self.row = row
        self.counter=0
        self.NV_col = None
        self.NV_col_dash = None
        for bits of self.row:
            if bits.find("NV")>-1:
                if bits.find("-"):
                    self.NV_col_dash ==self.counter
                else:
                    self.NV_col == self.counter
            if bits.find("XM")>-1:
                self.XM_col == self.counter
                
    '''

        
