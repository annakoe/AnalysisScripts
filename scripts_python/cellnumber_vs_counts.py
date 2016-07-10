### cellnumber_vs_counts.py v 1.0
### Anna Koeferle Nov 2015

###USAGE###
# python Cellnumber_vs_counts.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.colors as mplcolors
import matplotlib as mpl
import argparse
import sys

###Functions####

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

def replace_underscore_with_dash(dataset):
    dataset_renamed = dataset.rename(index=lambda x: re.sub(r"_", "-",x))
    return dataset_renamed

def guess_error_correction_method():
    if args.counts.split('/')[-1] == 'bayesian_corrected.csv' or args.counts.split('/')[-1] == 'bayesian_corrected_new.csv':
        method = 'Bayesian error correction'
    elif args.counts.split('/')[-1] == 'Dataframe_allsamples.txt':
        method = 'Barcodes grouped if < 4MM'
    elif args.counts.split('/')[-1] == 'Dataframe_allsamples_no_errors.txt':
        method = 'no error correction'
    else:
        method = 'unkown'
    return method

def restrict_to_samplesheet(mydataframe):
    samplefile = []
    new_df = pd.DataFrame()
    with open(args.samplesheet) as infile:
        for line in infile: #remove everything after the first underscroe and remove newline
            samplefile.append(re.sub(r"_.*", "", str(line.rstrip("\n"))))

    for sample in samplefile:
        if sample in mydataframe.index:
            print str(sample) + ' is in the dataset.'
            new_df = new_df.append(pd.Series(mydataframe.loc[sample, :]))  #append this index and all columns to dataframe
        else:
            print str(sample) + ' is NOT in the dataset.'

    return new_df


def plot_scatter(x_count, y_cellNumber, colorvar):
    fig,ax = plt.subplots() # use this to get ax object
    colors = ['midnightblue', 'darkorange']
    levels = [2, 3]
    cmap, norm = mplcolors.from_levels_and_colors(levels=levels, colors=colors, extend='max')

    myscatterplot = ax.scatter(x_count, y_cellNumber, c=colorvar, cmap=cmap, norm = norm, edgecolors='none', s=40) # plots the scatterplot

    colormap=myscatterplot.get_cmap()  #get the colormap used for the plot
    #use this to make proxy artists for the legend
    circles_artist=[mpl.lines.Line2D(range(1), range(1), color='w', marker='o', markersize=6, markeredgewidth=0, markerfacecolor=item) for item in colormap((np.array([2,3])-2)/1)]

    # mpl.lines.Line2D draws the circular object
    # using the colors from the list comprehension
    # colormap((np.array([2,3])-2)/1) # gets the colors out the colormap used for the plot. the array bit normalises the values to a range between 0-1 used for mapping

    ax.set_xlabel('Sum of gRNA counts', fontsize=16) # label for x axis
    ax.set_ylabel('Number of sorted cells', fontsize=16) # label for y axis

    method = guess_error_correction_method()
    title = fig.suptitle('All libraries - ' + str(method), fontsize=22) #figure title

    xdim = ax.get_xlim() # from the x object get the range of the data, returns a tuple
    xmin = xdim[0] # subset the tuple
    xmax = xdim[1]
    ydim = ax.get_ylim()
    ymin = ydim[0]
    ymax = ydim[1]
    identity_line = np.linspace(max(xmin, ymin),min(xmax, ymax)) #derive a  numpy array and get min and max comparing the values
    identity_plot = ax.plot(identity_line, identity_line, color="dimgrey", linewidth=2.0, label='x = y') # plot the identity line x=y

    line_artist = mpl.lines.Line2D([],[], color='dimgrey', linewidth=2.0, label='x = y') # draws the grey line in the legend
    ax.axis([xmin,xmax,ymin,ymax]) #rescale the axes so the identity line goes to edge of box

    line_legend = ax.legend(handles =[line_artist], loc = 4, frameon=False) #first legend is drawn inside the plot reads x= y
    # Add the legend for the line manually to the current Axes. This is the trick allowing us to plot two different legends!
    ax = plt.gca().add_artist(line_legend)
    # Create another legend for coloured points
    #lt.legend(handles=[line2], loc=4)

    circle_legend = plt.legend(circles_artist, ['2 cycles', '3 cycles'], loc = "center left", bbox_to_anchor = (1, 0.5), frameon =False, numpoints =1)
    plt.savefig('Sorted_cells_vs_gRNA_counts_' + str(method) + '.png', bbox_extra_artists=(circle_legend, title), bbox_inches='tight', format='png', dpi=300) ## bbox extra artist and bbox inches makes sure legend is not cut off figure
    plt.close()


###User input###

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--counts", help="dataframe containing counts per gRNA")
parser.add_argument("-n", "--cellnumber", help="dataframe containing counts of sorted cells")
parser.add_argument("-s", "--samplesheet", help="samplesheet containing one sample name per line, exactly as in dataframes")
args = parser.parse_args()

if not args.counts or not args.cellnumber:
   sys.exit("Missing inputfile!")

print "Reading data...."

delimiters = ['\t',', ', ',']

df_counts = find_delimiter(args.counts, delimiters)
if str(type(df_counts)) == "<type 'NoneType'>":
    print 'Delimiter not found in' + str(dfname_counts) + 'Please try again. Acceptable inputdataframes are tab-separated or comma-separated files.'
    sys.exit()

df_counts = remove_duplicates_from_df(remove_underscore_column_names(df_counts))
df_total_counts = pd.DataFrame(df_counts.sum(axis=0), columns = ['TotalCount'])
df_total_counts.index.name = None

number_of_cells = find_delimiter(args.cellnumber, delimiters)
if str(type(df_counts)) == "<type 'NoneType'>":
    print 'Delimiter not found in' + str(dfname_counts) + 'Please try again. Acceptable inputdataframes are tab-separated or comma-separated files.'
    sys.exit()

cell_count = remove_duplicates_from_df(remove_underscore_column_names(number_of_cells))
cell_count.drop(cell_count.columns[2:], axis=1, inplace=True)
cell_count = replace_underscore_with_dash(cell_count)
cell_count.index.name = None

frames = [df_total_counts, cell_count]
df_counts_cellnumber = pd.concat(frames, axis =1)
df_counts_cellnumber = df_counts_cellnumber.dropna(axis = 'index', how = 'any') #remove any rows that contain an NA
df_counts_cellnumber = restrict_to_samplesheet(df_counts_cellnumber) #makes sure only samples that are in the samplesheet are plotted in next step

plot_scatter(df_counts_cellnumber['TotalCount'], df_counts_cellnumber['CellCount'], df_counts_cellnumber['Cycles'])

print "Done!"
