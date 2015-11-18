#!/usr/bin/python
import os
import sys
import pdf2txt

class PDF_dump:
    '''
    REQUIRES PDFMINER (usuke Shinyama <yusuke at cs dot nyu dot edu>)
    http://www.unixuser.org/~euske/python/pdfminer/index.html#source
    create a pdf extractor object format = xml,text,tag,html
    >>> pdf = PDF_dump("paper2.pdf")
    >>> pdf.dump_file("html","outfile.txt")
    >>> txt = pdf.dump_obj("xml")
    >>>print txt
    '''
    def __init__(self,pdf_filename):
        self.filename = pdf_filename
        self.location = os.getcwd()

    def dump_file(self,tag="text",outfile=None):
        '''tag = xml, tag, text,html'''
        if outfile:
            outfil=self.location+"/"+outfile
        else:
            outfile=self.location+"/tmp.txt"
        flag = tag
        #path = self.location+"/"+self.filename
        cmd =  ["spacer","-t",flag,"-o",outfile, self.filename]
        try:
            pdf2txt.main(cmd)
        except IOError:
            print "Missing file!"

    def dump_obj(self,tag):
        '''returns a string object of pdf file in specified format'''
        self.dump_file(tag)
        tmpfile = self.location+"/tmp.txt"
        txtobj = open(tmpfile).read()
        os.system("rm "+tmpfile)
        return txtobj
        
    

                 


#os.system("python pdf2txt.py -t tag -o output4.txt gc_paper.pdf")
    #main(["pdf2txt.py","-t","tag","-o","output6.txt","paper2.pdf"])
