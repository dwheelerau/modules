#!/usr/bin/python
from linux_control import linux_controller
#wrapper for samtools to do things

sam_pth = "/home/dwheeler/software_tools/samtools-0.1.8/"

def sam2bam(sam_file,outfile):
    infile = sam_file
    #sam to sorted bam, memory eff, instrunctions from bowtie tutorial
    #http://bowtie-bio.sourceforge.net/tutorial.shtml
    #usage:sam2bam("VVadult_BS_m1_all_mapped.sam","testttt.bam")
    if outfile.find("."):
        #remove *.bam as this iwll be added automatically
        fix = outfile.split(".")
        outfile = fix[0]
    cmd1 = sam_pth+"./samtools view -bS -o tmp.bam " + infile
    print cmd1
    linux_controller(cmd1)
    print "sorting..."
    cmd2 = sam_pth+"./samtools sort tmp.bam " + outfile
    print cmd2
    linux_controller(cmd2)
    cmd3 = sam_pth+"./samtools index "+outfile+".bam"
    print cmd3
    print "indexing (creating *.bai file....."
    linux_controller(cmd3)
    print "check "+ outfile + ".bam and "+outfile+".bai" 

def sam2bam_alt(infile,outfile):
    #requires *fas file I think? this from google doc
    #convert to unsorted bam file
    cmd1= sam_pth+"./samtools import "+infile+ " tmp.bam"
    #sort this file
    cmd2 = sam_pth+"./samtools sort tmp.bam "+ outfile
    #index
    cmd3 = sam_pth+"./samtools index "+outfile
    linux_controller(cmd1)
    print cmd1
    linux_controller(cmd2)
    print cmd2
    linux_controller(cmd3)
    print cmd3
