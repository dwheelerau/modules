#a group of functions for reading and counting data in files
'''contains "strip_it" that strips chars from a file, "count_it" that
counts strings in a file, "read_it" that reads a specified number of
lines to the screen'''

from blast_tools import strip_it

def count_it(count_what, filename):
    '''takes a string to find and a filename searches through line by line
    and if finds string counts it. Returns number of finds as a int. If
    mulitple finds per line only counts once.
    '''
    find_it = str(count_what)
    counter = 0
    file_handle = open(filename,'r')

    for line in file_handle:
        if line.find(find_it)>=0:
            counter+=1
    file_handle.close()
    return counter

def read_it(number_of_lines,filename,start_from=1):
    '''Prints the number of lines asked for in the call from a file.
    flags = [-1] from the end, default is from the start [1]
    ****'all' prints the whole file to screen.****'''
    target = number_of_lines
    file_handle=open(filename,'r')
    counter = 0

    #count from begining of file
    if start_from == 1: #default
        for lines in file_handle:
            if counter<target:
                print lines

            counter+=1
    file_handle.close()
    
    #count from end of file
    if start_from == -1: #catch flag
        for lines in file_handle:
            counter+=1
        file_handle.close
        #reopen file
        file_handle = open(filename,'r')
        counter2=0
        for lines in file_handle:
            if counter2>counter-target:
                print lines
            counter2+=1

def count_lines(filename):
    '''simply counts the number of lines in a file'''
    file_handle = open(filename,"r")
    counter = 0
    for line in file_handle:
        counter+=1
    file_handle.close()
    return counter
        
        
