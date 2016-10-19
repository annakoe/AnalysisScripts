import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import collections as coll
import re
import pylab
import argparse
import sys

###Functions###

def remove_underscore_column_names(dataset):
    dataset_renamed = dataset.rename(columns=lambda x: re.sub(r"_.*", "",x))
    return dataset_renamed

def generate_samplefile(samplefilename):
    samplefile =[]
    with open(samplefilename) as infile:
        for line in infile: #remove everything after the first underscoe and remove newline
            samplefile.append(re.sub(r"_.*", "", str(line.rstrip("\n"))))
    return samplefile

def remove_3quarters_NAs(dataframe):
    df_with_NAs = dataframe.replace('0', np.nan)
    thresh = int(len(df_with_NAs.columns)/4*3)
    NA_removed =df_with_NAs.dropna(axis='index', thresh=thresh)
    return NA_removed

 def save_df_for_samplefile(samplefile, samplefilename): #this function calls all other functions
    df_samplefile = df_bayes[samplefile]
    df_samplefile = remove_3quarters_NAs(df_samplefile)
    df_samplefile = df_samplefile.fillna(0)
    df_samplefile.to_csv(str(samplefilename)+'_counts.csv')

###Script###
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdataframe", help="dataframe containing counts per gRNA after error correction")
parser.add_argument("-s", "--samplesheet", help="samplesheet containing one sample name per line, spelled exactly as in dataframes")
args = parser.parse_args()

if not args.inputdataframe or not args.samplesheet:
   sys.exit("Missing inputfile!")

print "Reading data...."


df_bayes = pd.DataFrame.from_csv(args.inputdataframe, header=0, sep=',', index_col=0)

samplefilename = args.samplesheet
samplefile =generate_samplefile(samplefilename)

save_df_for_samplefile(samplefile, samplefilename)
