#!/usr/bin/python
#FUCK!
from datetime import datetime

def date_file(leading_string):
    #return string with user supplied leader string wtih todays date
    #usage: date_file("bombus_seq")
    #returns leader with todays date: bombus_seq_2012-01-20
    currdate = datetime.today().isoformat()
    file_name = leading_string+"_"+currdate.split("T")[0]
    return file_name

#a = date_file("test")
    
    
