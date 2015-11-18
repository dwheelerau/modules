#!/usr/bin/python
'''general table tools, updated a bit better use of class structure
need to add support for using OGS file'''
class Row:
    '''defines a row with ability to split a line into key and val parts'''
    def __init__(self,key_index,row):
        self.row=row
        self.key_index = key_index

    def key_val(self):
        #returns key val pair as a string and string of \t delineated vals
        column=Column(self.row)
        key,val = column.split_col(self.key_index)
        return key,val
        
        
class Column:
    '''handles columns by split by tabs, returns split row'''
    def __init__(self,column):
        self.column = column.strip().split("\t")

    def split_col(self,split_at):
        #pops key out of column list and returns both as strings
        key = self.column.pop(split_at) #removes key from list!
        #col ="\t".join(self.column)
        return key,"\t".join(self.column).strip() #returns strings

class Table:
    def __init__(self,key,filename):
        self.filename = filename
        self.key = key
        self.tab_dic = {}

    def make_table_dic(self):
        #function for making a table dic with a specified key, removes key!
        for line in open(self.filename):
            row = Row(self.key,line)
            key_name,val_name = row.key_val()
            self.tab_dic[key_name]=val_name
        return self.tab_dic

class Table_joiner:
    '''a table joining class, collects info on filenames and keys, call
    self.join() to make a joined table. Access attribute self.missing_data list
    to get access to missing data (not found in one of second table but in
    the first can call directly Table_joiner(0,file1,0,file2,"outfile").join()
    Allows key reassignment if you want to make a new joined table using diff keys'''
    def __init__(self,key1,filename1,key2,filename2,outfile):
        self.key1 = key1
        self.key2=key2
        self.filename1=filename1
        self.filename2=filename2
        self.outfile = outfile
        self.missing_data = [] 

    def join(self):
        '''
        >>>Table_joiner(0,"infile.txt",0,"cpgoe.txt","outfile").join()
        >>>Table_joiner(0,"outfile",0,"blast_info2.txt","outfile2").join()
        '''
        #method that does the joining
        self.missing_data=[]
        outfile_handle = open(self.outfile,"w")
        #make table instances local that way can us one tablejoiner
        #change keys for example and re-join by calling method
        table1 = Table(self.key1,self.filename1).make_table_dic()
        table2 = Table(self.key2,self.filename2).make_table_dic()
        
        #take a random value and work out how many tabs for except keyerror
        table2_parts = len(table2[table2.keys().pop()].split("\t"))
        for key in table1:
            try:
                data=key+"\t"+table1[key]+"\t"+table2[key]+"\n"
            except KeyError:
                data=key+"\t"+table1[key]+(table2_parts*"\tNA")+"\n"
                self.missing_data.append(key)
            outfile_handle.write(data)
        if len(self.missing_data)!=0:
            open("check_me.txt","w").write("\n".join(self.missing_data))
            print "check file contains %s names"%len(self.missing_data)
        outfile_handle.close()

class Key_swap:
    '''a class that uses the ogs key file to take a old key and make
    replace this key with a new key, thus allowing it to be joined to
    another table by a common key'''
    
    def __init__(self,ogs_key_file,ogs_old_key,ogs_new_key):
        self.ogs_key_file = ogs_key_file
        #current key in the table that will be replaced
        self.ogs_old_key = ogs_old_key
        #new key will be this one
        self.ogs_new_key = ogs_new_key
        #self.translation_dic = {}

    def make_swap_dic(self):
        #make a dic with old key val as key and new key val as val
        translation_dic = {}
        for line in open(self.ogs_key_file):
            bits = line.strip().split("\t")
            translation_dic[bits[self.ogs_old_key]]=bits[self.ogs_new_key]
        return translation_dic

    def key_swap(self,tab_file,translation_dic,outfile):
        #swap first name in column(key) for another key from the ogs key file
        outfile_handle = open(outfile,"w")
        for line in open(tab_file):
            bits = line.strip().split("\t")
            old_key = bits[0]
            try:
                new_key = translation_dic[old_key]
                outfile_handle.write(new_key+"\t"+"\t".join(bits[1:])+"\n")

            except KeyError:
                pass
        outfile_handle.close()
            
        
if __name__=="__main__":
    #probably wont work!
    a=Table(0,"test.txt").make_table_dic()
    b=Table(0,"test2.txt").make_table_dic()
    Table_joiner(0,"test.txt",0,"test1.txt","outfile").join()
    out.close()
        
    
    
