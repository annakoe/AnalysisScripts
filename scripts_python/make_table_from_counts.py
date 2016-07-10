# makeTable_from_counts.py v 1.0
# Anna Koeferle OCt 2015, UCL
# This script takes the output of collapse_barcodes.py (either raw or no_orphans table) and formats the output into a dataframe for use with PCA and mageck scripts

import numpy as np
import pandas as pd
import fileinput
import sys
from time import time, clock
import pybedtools
from pybedtools import BedTool

####Global variables#####

### load library file into pybedtools BedTool object
### if using a library other than EMT5000 library, need to include path to library file here. Expects tab separated bed file of genomic target regions with associated target gene name
### e.g. chr "\t" start "\t" stop "\t" target gene name

EMT_lib = pybedtools.BedTool('/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts/EMT5000_library_regions.bed')

####Functions####

def load_file_to_dict(samplename):
    my_dict = {}
    with open(samplename, 'r') as infile:
        name = samplename
        for line in infile:
            lsplit = line.strip().split(', ')
            gRNA = lsplit[0]
            count = lsplit[1]
            my_dict[gRNA] = count
    return (my_dict, name) # tuple

# load_file_to_dict takes a filename as input and loads the corresponding file into a dictionary with keys: gRNA and values : gRNA count
# the inputfile is a comma-separated file of the form gRNA (chr:start-stop), count. This file is output by collapse_barcodes.py

def make_data_frame_from_dict(tuple_dict_name): #takes tuple
    my_dict = tuple_dict_name[0]
    samplename = tuple_dict_name[1]
    my_df = pd.DataFrame.from_dict(my_dict, orient='index')
    my_df.columns = [samplename.split("/")[-1].split("_")[0]]
    return my_df

# make_data_frame_from_dict takes a dictionary as input and makes it into a pandas dataframe using the keys as row indices

def Find_target_gene(my_dataframe, gRNA_library_with_genes):
    ###get index out of dataframe and into bed format
    gRNA_list =list(my_dataframe.index.values)

    gRNA_list_formatted = []

    for gRNA in gRNA_list:
        lsplit = gRNA.split (':')
        chrom = lsplit[0]
        startstop = lsplit[1].split('-')
        start = startstop[0]
        stop = startstop[1]
        gRNA_list_formatted.append( chrom + " " + start + " " + stop + " \n")


    gRNA_string = ' '.join(gRNA_list_formatted)
    gRNA_bed = BedTool(gRNA_string, from_string=True)
    gRNA_with_gene = gRNA_bed.intersect(EMT_lib, f=1, wa=True, wb=True)

    target_gene_name = [(f[6]) for f in gRNA_with_gene]
    my_dataframe.insert(0, "gene", target_gene_name)
    return my_dataframe


def catenate_dataframes():
    result = pd.DataFrame()  # make an empty pandas data frame

    for sample in sys.argv[1:]:   # loop over each file in the list supplied by the user
        new_df = make_data_frame_from_dict(load_file_to_dict(sample))  # get the dataframe returned when running the functions load_file_to_dict and make_data_frame_from_dict
        result = pd.concat([result, new_df], axis=1)  # update the empty data frame with the additional sample as a separate column

    result = Find_target_gene(result, EMT_lib) # inserts a column of target gene names for the gRNAs in the index column
    result.index.name = 'sgRNA'
    return result

###Script####

if (len(sys.argv) < 2):
    print "Missing inputfiles! Usage: python MakeTable_from_counts.py INPUTFILE1 INPUTFILE2 ... INPUTFILEn"
    sys.exit()

else:
    print "\nReading inputfile(s) into dictionary..."
    t0 = time() # Begin timer

    dataframe_output = catenate_dataframes()

    print "\nGeneratig outputfile..."
    dataframe_output.to_csv('Dataframe.txt', sep='\t', na_rep= '0') # na_rep tells what to output instead of NAs need to check what mageck can deal with

    timed = time() - t0
    print '\t Finished! Took', round(timed,3), 'seconds'
