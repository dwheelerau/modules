#!/usr/bin/python
#####################################################
#   example.py - a program to ....                  #
#                                                   #
# Author: Dave Wheeler                              #
#                                                   #
# Purpose:extract species id from ncbi blastx       #
# draw a histogram of top 10 species                #
# Usage:                                            #
#####################################################
import sys
from pandas import Series, DataFrame
import os
import matplotlib.pyplot as plt

def get_genus(hit_line):
    """get genus hit info"""
    genus = extract_genus(hit_line)
    #else:
    #genus = "No hit"
    return genus

def get_species(hit_line):
    """split table by tab and get spec hit info"""
    species = extract_species(hit_line)
    return species

def extract_genus(info):
    """extract Drosophila from [Drosophila erecta] tags"""
    genus_spec = info[info.find('[')+1:info.find(']')]
    return genus_spec.split(' ')[0]

def extract_genus(info):
    """extract Drosophila from [Drosophila erecta] tags"""
    genus_spec = info[info.find('[')+1:info.find(']')]
    return genus_spec.split(' ')[1]

def get_hist_data(alist,limit=10):
    """make a freq series from list using pandas, limit to most freq 10"""
    genus_freq = Series(alist).value_counts()
    genus_freq.sort(ascending=False) #biggest to smallest
    return genus_freq[:limit]

if __name__ == "__main__":
    table = sys.argv[1]
    if len(sys.argv) >2:
        outfile = sys.argv[2]
    else:
        outfile = "test.png"
    infile = open(table)
    geni = []
    header = infile.next()
    for line in infile:
        genus = get_genus(line)
        geni.append(genus)
    infile.close()
    genus_freq = DataFrame(get_hist_data(geni,8))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    p = genus_freq.plot(kind="bar",legend=False,color='grey',rot=50,axes=ax)
    p.set_title("Best blastx hit by genus")
    p.set_ylabel("Frequency")
    #p.set_xlabel("Genus")
    figs = p.get_figure()
    figs.savefig("%s/%s"%(os.getcwd(),outfile))








