#plot_counts_vs_number_of_reads.py v 1.0.0
# Anna Koeferle, UCL, November 2015
#takes as input a dataframe holding number of reads per gRNA, a dataframe holding number of counts per gRNA and a samplesheet, listing the samples to be processed.

#####USAGE######

# python plot_counts_vs_number_of_reads.py [Dataframe-Counts.csv] [Dataframe-NumberOfReads.csv] [Path/to/list_of_Samplenames]


from __future__ import division
import numpy as np
import pandas as pd
import sys
import re
import matplotlib.pyplot as plt
import pylab as pylab


####Functions####

def find_delimiter(filename, delimiters):
    #checks if file is csv or tsv and loads file into pandas dataframe
    for delim in delimiters:
        df_counts = pd.DataFrame.from_csv(filename, header=0, sep= delim, index_col=0)
        if len(df_counts.columns) > 0:
            print 'Delimiter found!'
            break
        else:
            print 'Have not found the correct delimiter yet. Still searching....'
            df_counts = None
    return df_counts

def remove_underscore_column_names(dataset):
    dataset_renamed = dataset.rename(columns=lambda x: re.sub(r"_.*", "",x))
    return dataset_renamed

def remove_duplicates_from_df(dfname):
    new_df = dfname.T.groupby(level=0).first().T
    return new_df

def plot_scatterplot(sample):
    new_df = pd.concat([df_reads[sample], df_counts[sample]], axis =1, keys=['reads', 'counts'])
    fig = plt.figure()
    plt.scatter(new_df['reads'], new_df['counts'], c='steelblue')
    plt.xlabel('number of reads', fontsize=16)
    plt.ylabel('number of counts (unique molecules)', fontsize=16)
    fig.suptitle(str(sample), fontsize=24)
    plt.savefig(str(sample) + '.png', format='png', dpi=300)
    plt.close()

#####Checking user input#####

if (len(sys.argv) <> 4):
    print "Missing inputfile(s)! Usage: python plot_counts_vs_number_of_reads.py [Dataframe-Counts] [Dataframe-NumberOfReads] [Samplesheet]"
    sys.exit()

else:
    dfname_counts = str(sys.argv[1])
    dfname_reads = str(sys.argv[2])
    name_samplefile = str(sys.argv[3])

    delimiters = ['\t',', ', ',']

    ####load the count dataframe:
    df_counts = find_delimiter(dfname_counts, delimiters)
    if str(type(df_counts)) == "<type 'NoneType'>":
        print 'Delimiter not found in' + str(dfname_counts) + 'Please try again. Acceptable inputdataframes are tab-separated or comma-separated files.'
        sys.exit()

    df_counts = remove_duplicates_from_df(remove_underscore_column_names(df_counts))

    ###load the read dataframe:
    df_reads = find_delimiter(dfname_reads, delimiters)
    if str(type(df_reads)) == "<type 'NoneType'>":
        print 'Delimiter not found in' + str(dfname_reads) + 'Please try again. Acceptable inputdataframes are tab-separated or comma-separated files.'
        sys.exit()

    df_reads = remove_duplicates_from_df(remove_underscore_column_names(df_reads))

    ####load the samplesheet into a list:
    samplefile =[]
    with open(name_samplefile) as infile:
        for line in infile: #remove everything after the first underscoe and remove newline
            samplefile.append(re.sub(r"_.*", "", str(line.rstrip("\n"))))

    for sample in samplefile:
        try:
            plot_scatterplot(sample)
        except Exception, msg:
            print str(msg)
