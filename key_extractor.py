#!/usr/bin/python
import os

class Key_extractor:
    def __init__(self,target, file_name="OGS_v1.2_official_id_map.txt",OS="linux"):
        if OS=="linux":
            CONNECT = "/"
        else:
            CONNECT = "\\"
        self.location = os.getcwd()
        #self.file_name = self.location+CONNECT+file_name
        #self.file_name = "/home/dwheeler/scripts/OGS_v1.2_official_id_map.txt"
        self.file_name = file_name
        self.target = target
        self.data = self.get_info()
                
    def get_info(self):
        handle = open(self.file_name)
        data=None
        for line in handle:
            if line.find(self.target)!=-1:
                data=line.split("\t")
                return data
        if data==None:
            return ["not found","not found","not found","not found",\
                        "not found","not found"]

    
if __name__=="__main__":
    test = Key_extractor("NV22491").get_info()
    print test
