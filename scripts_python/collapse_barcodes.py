# collapse_barcodes.py v 1.1
# Anna Koeferle Dec 2015, UCL
# adapted from CollapseTCRs.py v1.2 by James Heather and Katharine Best

### BACKGROUND ###
# Makes use of random molecular barcode sequences to error-correct high-throughput sequencing data.
# The barcodes are random nucleotides (N) added to the library amplicon prior to amplification.
# Instead of counting reads/molecules (which is problematic due to PCR amplification bias), the co-occurence of different barcodes with a sequence of interst is counted (count barcodes instead of reads).
# This allows us to correct for PCR amplification bias and PCR error/sequencing error.
# The script does the following:
# 1. Import input file of the form (c1) identifier (c2) gRNA chr:start-end OR DNA sequence (c3) barcode (varying lengths)
# 2. The script then sorts according to (c2) gRNA.
# 3. Within each gRNA, the barcodes are grouped according to sequence similarity
# If a barcode is within x edit distance from the previous one group them together
# 4. barcodes in groups are re-written/corrected to the sequence of the most common member of the group
# 5. export a file that has the corrected barcode sequence in c4

### INPUT ###
# takes a tab-delimited file consisting of 3 columns:
# column1: Read identifier
# column2: gRNA from the library read was aligned back to, either as DNA sequence or coordinates of the form chr:start-end
# column3: barcode (various lengths N)
# Run: python CollapseBarcodess.py FILENAME.n12

### OUTPUT ###
### Outputs a gRNA frequency file = gRNA and its counts
## can output a .dict file = the collapsed dictionaries used to calculate the gRNA frequencies, by default commented out

### USAGE ###
# python CollapseBarcode.py INPUTFILENAME MAX_NUMBER_OF_ERRORS

#### PACKAGES ####

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

##### Setting functions and empty dictionaries for later use
def dodex():
    return coll.Counter()

def breakdown(etc):
  # Used to break a given dcr_etc (i.e. what is stored in a given line of the inputfile) into its components
  # Splits on '|', to avoid breaking up within identifiers or fastq quality strings

  return re.findall(r"[\w,\<\=\>\-\;\:\?\/\.\@\# ]+", str(etc))
  # breakdown[0] = gRNA / [1] = ID

def guide_frequency(collections_dictionary):
    output=[]
    for i in collections_dictionary:
        clustered = coll.Counter()
        for k,v in collections_dictionary[i].most_common():
            clustered[k] += v
        result = str(i) + ", " + str(len(clustered)) + "\n" # len of dictionary returns number of key:values stores output in string
        output.append(result)

    return output

def remove_orphan_barcodes(diction):
    #removes barcodes that are seen only once, even after error correction
    dict_morethan1 = coll.defaultdict(list)

    for x in diction: #loop over the gRNAs
        counter = coll.Counter() # set an empty counter
        for k,v in diction[x].most_common():  # loop over the barcodes and counts for each gRNA
            if v > 1:    # if the count is greater than 1
                counter[k] += v   # add the barcode and counter to the counter dictionary
        dict_morethan1[x] = counter   # append to the gRNAs

    return dict_morethan1


#### Importing data from input###

total_number_of_reads = 0 # count_input_lines

if (len(sys.argv) <> 2):
    print "Missing inputfile! Usage: python collapse_barcodes_editdist0.py INPUTFILENAME"
    sys.exit()
else:
    seqfilename = str(sys.argv[1])    # stores name of user-supplied file in variable seqfilename
    distance_threshold = int(sys.argv[2]) #stores user-supplied max number of edit distances for barcodes to be grouped together
    start_time = time() # set time at start

    print "\nReading", str(seqfilename), "into dictionary..."
    t0 = time() # Begin timer

### TAB- DELIMITED POSITIONS IN THE INPUTFILE ARE: ###

# Read identifier -[0]-
# gRNA -[1]-
# Barcode-[2]-

    dict_seqfile = coll.defaultdict(list)   #make an empty collections default dictionary

    with open(seqfilename) as infile:
        for line in infile:
            total_number_of_reads += 1
            lsplit = line.strip().split('\t')
            readID = lsplit[0]
            gRNA = lsplit[1]
            barcode = lsplit[2] #this includes barcodes of varying lengths

            values = gRNA + "|" + readID
            dict_seqfile[barcode].append(values) # this now is a dictionary with key=barcode and value=a list of gRNA and readIDs separated by "|" all with the same key, e.g. looks like 'TGTGGGCGACGGGG': ['chr10:31609045-31609068|>D00623:46:H7JVFBCXX:1:1101:5035:69776.', 'chr10:31609045-31609068|>D00623:46:H7JVFBCXX:1:1102:10781:28677.', 'chr10:31609045-31609068|>D00623:46:H7JVFBCXX:1:1102:9779:13200.', 'chr10:31609045-31609068|>D00623:46:H7JVFBCXX:1:1103:3007:69834.', 'chr10:31609045-31609068|>D00623:46:H7JVFBCXX:1:1104:1364:82596.', 'chr10:31609045-31609068|>D00623:46:H7JVFBCXX:1:1104:15350:25118.']

    timed = time() - t0
    print 'The total number of reads for this sample was:' +str(total_number_of_reads) +'.'

#### now take the inputfile and derive something that looks like dcr_collapsed dictionary of Jamie Heather:

print "Counting barcodes associated with each gRNA..."

dict_collapsed = coll.defaultdict(dodex)
t0 = time() # Begin timer


for key in dict_seqfile:  #loop throught the barcodes, i.e. the key. like saying 'for every barcode in the dictionary':
    if len(key) >= 14:   # exclude barcodes that are less than 14 bases in length
        for value in dict_seqfile[key]:          # loop through the values (gRNA | identifier)
            gRNA_seq = breakdown(value)[0] # prints only the value correspondign to gRNA
            dict_collapsed[gRNA_seq][key] += 1 # appends to dict_collapsed the gRNA sequence, the barcode = key and its count. In fact, does the counting

timed = time() - t0

#print'InputData collapsed into dictionary: '
#print dict_collapsed

print '\t Finished! Took', round(timed,3), 'seconds'

###### dict_collapsed looks like this with test datasset: {gRNA: Counter({'barcode key':count, 'barcode':count}),
#defaultdict(<function dodex at 0x103000de8>, {'chr10:31609167-31609190': Counter({'TTA': 1888, 'TTATTT': 985, 'AAGTAA': 752, 'CCCCCG': 8, 'CTCTTC': 4, 'CTTATG': 4,

print "Removing orphan barcodes..."

t0 = time() #begin timer

try:
    dict_no_orphans = remove_orphan_barcodes(dict_collapsed)
except Exception, msg:
    print str(msg)


timed = time() - t0
print '\t Finished! Took', round(timed,3), 'seconds'

###### calculate the frequencies of gRNAs based on length of the dictionaries
print "Computing gRNA frequencies..."

try:
    frequency_clustered = guide_frequency(dict_collapsed)
except Exception, msg:
    print str(msg)

with open(seqfilename.split(".")[0]+ "_frequency_raw", 'w') as outfile:
    for item in frequency_clustered:
        outfile.write(item)

########### gRNA counts for dict_threshold

try:
    frequency_no_orphans = guide_frequency(dict_no_orphans)  # run the function guide_frequency on the dictionaries generated above and store output in a string
except Exception, msg:
    print str(msg)


with open(seqfilename.split(".")[0] + "_frequency_no_orphans", 'w') as outfile: #open a csv file to write to
    for item in frequency_no_orphans: #loop through the elements of the list, and print elements of list to file
        outfile.write(item)

print "\tDONE!"

#print "Outputting dictionaries to file...."

#json.dump(dict_clustered, open(seqfilename.split(".")[0]+str("_collapsed_dict_readcountnorm"), 'w'))
#json.dump(dict_no_orphans, open(seqfilename.split(".")[0]+str("_collapsed_no_orphans_dict_readcountnorm"), 'w'))
#json.dump(dict_threshold, open(seqfilename.split(".")[0]+str("_collapsed_threshold_dict_readcountnorm"), 'w'))
