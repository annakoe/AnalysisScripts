### plot_LogFC_vs_LogFC.py v 1.0
### Anna Koeferle Dec 2015

###USAGE###
# python plot_LogFC_vs_LogFC.py [list of file to be compared]

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import collections as coll
import re
import pylab
import pybedtools
from matplotlib import gridspec
import sys

def read_filenames():
    list_of_filenames =[]
    for sample in sys.argv[1:]:   # loop over each file in the list supplied by the user
        filename = sample.split('/')[-1].split('_DESeq2')[0] #splits the filename at all '/' and takes the last element, then uses everything before '_DESeq2'
        list_of_filenames.append(filename)
    return list_of_filenames


def read_dataframes():
    result = pd.DataFrame()  # make an empty pandas data frame
    for samplefilepath in sys.argv[1:]:   # loop over each file in the list supplied by the user
        df_DESeq2 = pd.DataFrame.from_csv(samplefilepath, header=0, sep=',', index_col=0) #loads the DESEq2 output into dataframe
        result = pd.concat([result, df_DESeq2["log2FoldChange"]], axis=1)  # update the empty data frame with the column of Log2Fold Changes from the DESeq2 output
        result.columns = result.columns.str.replace('log2FoldChange', str(samplefilepath.split('/')[-1].split('_DESeq2')[0]))
    return result

def generate_outputfilename():
    string = ''
    for file in list_of_filenames: #loop over filenames and make a string of these
        string += file
    return string

def make_dataframe_for_plotting(x_L2FC, y_L2FC):
    series =[x_L2FC,y_L2FC] #put the two columns into a list
    my_df=pd.concat(series, axis=1) #concatenate the two columns in the list into a dataframe
    my_df.fillna(0) #replace all NAs by 0
    return my_df


def plot_the_grid(filenames, dataframe_L2FC):
    num_iter = len(filenames) #number of iterations of the loop
    dim = len(filenames)-1 #dimensions of the plotting grid

    gs = gridspec.GridSpec(dim, dim) #set a grid to the dimension of the matrix of subplots
    iterstart =0
    for x in range(0, num_iter): #loop over the x dimension (loop over columns)
        row_coord =0
        iterstart +=1
        if iterstart >1: #after the first loop over the first column is complete this will be executed
            ax.set_xlabel(filenames[x-1], fontsize=8) #add the xlabel from the previous iteration, this labels only the last plot in a row of plots
            ax.get_xaxis().set_tick_params(direction='out', labelbottom='on', labelsize=8)

        for y in reversed(range(iterstart, num_iter)): #loop over the rows, the start of the yrange is increased at every iteration meaning for every column we move to the right of the grid, there will be one fewer row plotted
            print 'Plotting ' + str(filenames[x]) + ' at pos ' + str(x) + ' versus ' + str(filenames[y]) + ' at pos ' + str(y)
            print 'Plot position: Row = ' + str(row_coord) + ' and Column = ' + str(x)

            #ax = plt.subplot(gs[row, column])
            ax = plt.subplot(gs[row_coord, x]) #where to place this subplot on the grid gs
            new_dataframe = make_dataframe_for_plotting(dataframe_L2FC[filenames[x]], dataframe_L2FC[filenames[y]]) #put the columns to plot against each other in a temporary dataframe

            #plot_scatter(ax, dataframe_L2FC[filenames[x]], dataframe_L2FC[filenames[y]], filenames[x], filenames[y], x, row_coord) #plot the subplot
            ax.scatter(new_dataframe.iloc[:,[0]], new_dataframe.iloc[:,[1]], c='black', edgecolor='none', s=8) # plots the scatterplot x against y
            ax.set_xlim([-10, 10])
            ax.set_ylim([-10, 10])

            xmin,xmax = ax.get_xlim()
            ymin,ymax = ax.get_ylim()
            ax.set_aspect(abs(xmax-xmin)/abs(ymax-ymin)) #make a square plot

            ax.spines['right'].set_visible(False) #remove the right axis
            ax.spines['top'].set_visible(False) #remove the top axis

            plt.axhline(0, color='black', ls='--', lw=1) #draws a horizontal black line through 0
            plt.axvline(0, color='black',ls='--', lw=1) #draws a vertical black line through 0

            if x == 0: #for subplots in the first row of the gid only
                ax.set_ylabel(filenames[y], fontsize=9) # label for y axis
                ax.get_yaxis().set_tick_params(direction='out', labelsize=8)
            else:
                ax.set_yticklabels([]) #no y labels for all other positions on the grid

            ax.yaxis.set_ticks_position('left')
            ax.get_yaxis().set_tick_params(direction='out')

            ax.get_xaxis().set_tick_params(direction='out', labelbottom='off', top='off') #no x axis ticklabels on any of the plots

            row_coord += 1 #increment the row coordinate by one

    plt.savefig('All_by_all_L2FC_' + generate_outputfilename() + '.png', bbox_inches='tight', format='png', dpi=300)
    plt.close()



if (len(sys.argv) < 2):
    print "Missing inputfiles! Usage: python plot_LogFC_vs_LogFC.py INPUTFILE1 INPUTFILE2 ... INPUTFILEn"
    sys.exit()

print "\nReading inputfile(s) into dataframe..."

list_of_filenames = read_filenames()
dataframe_L2FC = read_dataframes()
plot_the_grid(list_of_filenames, dataframe_L2FC)
