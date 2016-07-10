#make_csv_4Bayes.py v 1.0
#written by Anna Koeferle, 2015
#adapted in large parts from  the script CollapseTCRs.py by Jamie Heather (can be found at https://github.com/JamieHeather/tcr-analysis/blob/master/CollapseTCRs.py)

#This script takes in a tab-separated file containing  [0] read ID [1] gRNA chr:start-end [2] barcode 5primer and 3 prime fused as input
#This script ouptuts a csv file to input into the bayesian PCR error correction script written by James E. Barrett
#The output is a csv file with columns:  [0] gRNA chr:start-end,  [1] barcode (14 bp 5' and 3' barcode combined, [2] read count for this barcode gRNA combination (no error correction)

from __future__ import division
import collections as coll
import sys
import Levenshtein as lev
import re
from operator import itemgetter
from time import time, clock
import json
import signal
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO
import numpy as np
import matplotlib.pyplot as plt
import pylab as pylab


def dodex():
    return coll.Counter()

def breakdown(etc):
  # Splits on '|', to avoid breaking up within identifiers or fastq quality strings

  return re.findall(r"[\w,\<\=\>\-\;\:\?\/\.\@\# ]+", str(etc))
  # breakdown[0] = gRNA / [1] = ID

def read_data(seqfilename):
    dict_seqfile = coll.defaultdict(list)   #make an empty collections default dictionary

    with open(seqfilename) as infile:
        for line in infile:
            lsplit = line.strip().split('\t')
            readID = lsplit[0]
            gRNA = lsplit[1]
            barcode = lsplit[2]
            values = gRNA + "|" + readID
            dict_seqfile[barcode].append(values)

    return dict_seqfile

def count_barcodes(dict_seqfile):

    dict_collapsed = coll.defaultdict(dodex)

    for key in dict_seqfile:  #loop throught the barcodes, i.e. the key. like saying 'for every barcode in the dictionary':
        if len(key) >= 14:   # exclude barcodes that are less than 14 bases in length
            for value in dict_seqfile[key]:          # loop through the values (gRNA | identifier)
                gRNA_seq = breakdown(value)[0] # prints only the value correspondign to gRNA
                dict_collapsed[gRNA_seq][key] += 1 # appends to dict_collapsed the gRNA sequence, the barcode = key and its count. In fact, does the counting

    return dict_collapsed

#write the dictionary into a csv file

def write_dict_to_csv(dict_clustered_name, outfilename):
    outfile = open( outfilename.split(".")[0]+str("_no_error_correction.csv"), 'w' )

    for key, value in dict_clustered_name.items():  # iterate through the coll.dictionary
        for item in value.iteritems(): #iterate through the counter
            outfile.write( str(key) + "," + ','.join(map(str, item)) + '\n')

if (len(sys.argv) <> 2):
    print "Missing inputfile! Usage: python make_csv_james_model.py INPUTFILENAME"
    sys.exit()
else:
    seqfilename = str(sys.argv[1])    # stores name of user-supplied file in variable seqfilename
    start_time = time() # set time at start

    print "\nReading", str(seqfilename), "into dictionary..."


dict_seqfile = read_data(seqfilename)
dict_collapsed = count_barcodes(dict_seqfile)

print '\t Finished! Writing to csv'

write_dict_to_csv(dict_collapsed, seqfilename)

print '\t DONE!'
