#!/usr/bin/python
#file utilities to read and write files in my share data format
####################### like this ########## 
#Owner = Dave
#AMP.AX	988	5.37	5305.56
#ARG.AX	1344	6.07	8158.08
#Owner = Us
#AIO.AX	660	1.60	1056.00
#BHP.AX	80	45.03	3602.40
#end
#####################

from Share import *

class Share_file:
    '''a class for reading and writing my standard share files
    provides access to data through self.owner_shares = owner name
    as key and a list of Share instances. Also provides access to
    self.owner_value dict which contains owner =key and holding old value
    for comparison with new value as a Decimal object.
    The new value can be obtained by calling
    Decimal(share.get_value()*share.number) after iterating through the
    list using ie "share" -> share.value, share.ticker etc
    '''
    
    def __init__(self,filename):
        self.filename = filename
        #could pop these here by calling read_file? Not sure
        self.owner_shares={}
        self.owner_value={}
        
    def read_file(self):
        '''reads standard share file and returns a dictionary of owners
        with shareobjects and also returns a dictionary of old total values
        from the file, this could be used to check for increases or dec
        '''
        file_handle = open(self.filename)
        #owner_shares={}
        #owner_value={}
        for line in file_handle:
            line = line.strip()
            if line.find("Owner")>-1:
                #share_list = []
                owner = line.split("=")[1].strip()
                #setup_list
                self.owner_shares[owner]=[]
            elif line.find("#end")!=-1:
                break
                #DOES not return anything, needs to be accessed throguh instance?
                #not sure if this is desirable or not, probably good?
                #return self.owner_shares, self.owner_value
            else:
                ticker,number,purchase_price,old_val = line.split("\t")
                self.owner_shares[owner].append(Share(ticker,number,purchase_price))
                self.owner_value[owner]=Decimal(old_val)
        file_handle.close()
        
    def write_file(self,alt_filename=None):
        # write out on the proper format this is why set up dics as self!
        assert len(self.owner_shares)!=0
        backup_handle=open(self.filename+".bak","w")
        file_handle = open(self.filename).read()
        backup_handle.write(file_handle)
        backup_handle.close()
        if alt_filename==None:
            outfile_handle = open(self.filename,"w")
        else:
            outfile_handle=open(alt_filename,"w")
        for name in self.owner_shares:
            #returns share holder name, use that to access share list
            outfile_handle.write("Owner = "+name+"\n")
            for share in self.owner_shares[name]:
                #could just call directly when unpacking string formating?
                ticker = share.ticker
                number = share.number
                price = share.value
                value = Decimal(share.get_value()*number)
                #write out in correct format
                outfile_handle.write("%s\t%s\t%s\t%s\n"%(ticker,number,price,value))
        #close the file with a end of file tag
        outfile_handle.write("#end")
        outfile_handle.close()
                        
if __name__ == "__main__":
    test= Share_file("/home/dwheeler/Desktop/test.txt")
                                
