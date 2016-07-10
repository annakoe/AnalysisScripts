#plot_counts_vs_counts.py v 1.0.0
# Anna Koeferle, UCL, November 2015
#takes as input two dataframes (column: samplename, rows: gRNA counts, index: gRNA as chr:start-stop) holding number of counts per gRNA, and a samplesheet, listing the samples to be processed


from __future__ import division
import numpy as np
import pandas as pd
import sys
import re
import matplotlib.pyplot as plt
import pylab as pylab

#####USAGE######

# python plot_counts_vs_counts_totalcountnorm.py [Dataframe-Counts-xaxis.csv] [Dataframe-Counts-yaxis.csv] [Path/to/list_of_Samplenames]

####Functions####

def find_delimiter(filename, delimiters):
    #checks if file is csv or tsv and loads file into pandas dataframe
    for delim in delimiters:
        df_counts = pd.DataFrame.from_csv(filename, header=0, sep= delim, index_col=0) #read in file using first delimiter in list of possible deliminters
        if len(df_counts.columns) > 0: #check if using this delimiter has split the file correctly into columns
            print 'Delimiter found!'
            break
        else:
            print 'Have not found the correct delimiter yet. Still searching....'
            df_counts = None #set the delimiter to None if it was incorrect
    return df_counts

def remove_underscore_column_names(dataset):
    dataset_renamed = dataset.rename(columns=lambda x: re.sub(r"_.*", "",x)) #removes everything in sample name after first occurence of an underscore
    return dataset_renamed

def remove_duplicates_from_df(dfname):
    new_df = dfname.T.groupby(level=0).first().T #removes duplicate colunns from df
    return new_df

def plot_scatterplot(sample):
    new_df = pd.concat([df_counts1[sample], df_counts2[sample]], axis =1, keys=['counts1', 'counts2'])
    fig,ax = plt.subplots() # use this to get ax object
    ax.scatter(new_df['counts1'], new_df['counts2'], c='midnightblue', label = None) # plots the scatterplot
    if dfname_counts1.split('/')[-1] == 'bayesian_corrected.csv':
        ax.set_xlabel('counts (Bayesian error correction)', fontsize=16) # label for x axis
    elif dfname_counts1.split('/')[-1] == 'Dataframe_allsamples.txt':
        ax.set_xlabel('counts (PCR error correction max 4mm)', fontsize=16) # label for x axis
    elif dfname_counts1.split('/')[-1] == 'Dataframe_allsamples_no_errors.txt':
        ax.set_xlabel('counts (no error correction)', fontsize=16) # label for x axis
    else:
        ax.set_xlabel('Unknown', fontsize=16) # label for x axis

    if dfname_counts2.split('/')[-1] == 'bayesian_corrected.csv':
        ax.set_ylabel('counts (Bayesian error correction)', fontsize=16) # label for y axis
    elif dfname_counts2.split('/')[-1] == 'Dataframe_allsamples.txt':
        ax.set_ylabel('counts (PCR error correction max 4mm)', fontsize=16) # label for y axis
    elif dfname_counts2.split('/')[-1] == 'Dataframe_allsamples_no_errors.txt':
        ax.set_ylabel('counts (no error correction)', fontsize=16) # label for y axis
    else:
        ax.set_ylabel('Unknown', fontsize=16) # label for y axis

    fig.suptitle(str(sample), fontsize=24) #figure title
    xdim = ax.get_xlim() # from the x object get the range of the data, returns a tuple
    xmin = xdim[0] # subset the tuple
    xmax = xdim[1]
    ydim = ax.get_ylim()
    ymin = ydim[0]
    ymax = ydim[1]
    identity_line = np.linspace(max(xmin, ymin),min(xmax, ymax)) #derive a  numpy array and get min and max comparing the values
    ax.plot(identity_line, identity_line, color="darkorange", linewidth=2.0, label='x = y') # plot the identity line x=y
    ax.axis([xmin,xmax,ymin,ymax]) #rescale the axes so the identity line goes to edge of box
    ax.legend(loc =4, frameon =False)
    plt.savefig(str(sample) + '.png', dpi=300)
    plt.close()

#####Checking user input#####

if (len(sys.argv) <> 4):
    print "Missing inputfile(s)! Usage: python plot_counts_vs_counts_totalcountnorm.py [Dataframe1-Counts] [Dataframe2-Counts] [Samplesheet]"
    sys.exit()

else:
    dfname_counts1 = str(sys.argv[1])
    dfname_counts2 = str(sys.argv[2])
    name_samplefile = str(sys.argv[3])

    delimiters = ['\t',', ', ',']

    ####load the count dataframe:
    df_counts1 = find_delimiter(dfname_counts1, delimiters)
    if str(type(df_counts1)) == "<type 'NoneType'>":
        print 'Delimiter not found in' + str(dfname_counts1) + 'Please try again. Acceptable inputdataframes are tab-separated or comma-separated files.'
        sys.exit()

    df_counts1 = remove_duplicates_from_df(remove_underscore_column_names(df_counts1))

    ###load the read dataframe:
    df_counts2 = find_delimiter(dfname_counts2, delimiters)
    if str(type(df_counts2)) == "<type 'NoneType'>":
        print 'Delimiter not found in' + str(dfname_counts2) + 'Please try again. Acceptable inputdataframes are tab-separated or comma-separated files.'
        sys.exit()

    df_counts2 = remove_duplicates_from_df(remove_underscore_column_names(df_counts2))

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
