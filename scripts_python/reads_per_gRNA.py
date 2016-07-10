# Reads_per_gRNA.py v 1.0
# Anna Koeferle Nov 2015, UCL
# Reads_per_gRNA.py counts the number of reads associated with each gRNA.
# This script takes the a tab-separated file as input:
# the Inputfile has the following columns: [0] read ID [1] gRNA chr:start-stop [2] 14 nt barcode
# the output is a csv file: [0] gRNA chr:start-stop , [1] number of reads

from __future__ import division
import collections as coll
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import pylab as pylab

####### Functions ######

def read_data(seqfilename):
    dict_seqfile = coll.defaultdict(list)   #make an empty collections default dictionary

    with open(seqfilename) as infile:
        for line in infile:
            lsplit = line.strip().split('\t')
            readID = lsplit[0]
            gRNA = lsplit[1]
            barcode = lsplit[2]
            dict_seqfile[gRNA].append(readID)
    return dict_seqfile


def write_dict_to_csv(dict_name, outfilename):
    with open(outfilename, 'w') as outfile:
        for key in dict_name:
            outfile.write(str(key) + ', ' + str(len(dict_name[key])) + '\n')


def sanity_check(dict_seqfile):
    sanity = {}  ### calculate total number of reads in sample and crosscheck
    for key in dict_seqfile:
        sanity[key] = len(dict_seqfile[key])

    sumcheck = 0
    for key in sanity:
        sumcheck += sanity[key]

    return sumcheck

######Script#####

if (len(sys.argv) <> 2):
    print "Missing inputfile! Usage: python Reads_per_gRNA.py INPUTFILENAME"
    sys.exit()
else:
    seqfilename = str(sys.argv[1])
    outfilename = seqfilename.split(".")[0]+str("_read_count")

    print "\nReading", str(seqfilename), "into dictionary...and printing to csv..."

    dict_seqfile = read_data(seqfilename)

    write_dict_to_csv(dict_seqfile, outfilename)

    total_number_of_reads = sanity_check(dict_seqfile)

    print "Done! The total number of reads for this sample were: " + str(total_number_of_reads)
