#!/usr/bin/python
#import os this is outof date
import subprocess
from os import listdir
from os.path import isfile, join
from os import getcwd
'''
Some functions for controlling the terminal example:

>>>run_ls("test","*.txt")
>>> for names in open("test"):
	names = names.strip() #remove trailing new-line
	linux_controller("cp "+names+" tst")#copies files to a directory

###Update to subprocess###
subprocess.check_output("ls test*; exit 0",shell=True,stderr=subprocess.STDOUT)

'''

def run_ls(file_name,pattern=""):
    '''sends ls to file_name, opt pattern ie ls pattern ie[*.doc]
    usage: run_ls("file_list","*.txt")'''
    #os.system("ls " + pattern + "> "+file_name)
    cmd = 'ls '+ pattern + '> ' + file_name +' ;exit 0'
    subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
    
def linux_controller(command_string):
    '''generic wrapper for terminal conmands, takes command string
    >>> linux_controller("echo 'test' > test") prints "test" to file
    >>>for names in file_list:
    >>>     linux_controller("'cp '+names
    Returns 1 for failed
    '''
    #os.system(command_string)
    a= subprocess.call(command_string,shell=True,stderr=subprocess.STDOUT)
    return a


def file_list(directory, pattern=""):
    '''returns a list of filenames based on a pattern'''
    onlyfiles = [
        f for f in listdir(directory)
        if isfile(join(directory, f))
        if f.find(pattern) > -1]
    return onlyfiles



#test
#run_ls("test","test*")
#a = linux_controller("touch tttt")
#b = file_list("test*")

