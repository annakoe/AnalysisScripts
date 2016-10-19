### Plot_Log2FC_along_chr.py v 1.0
### Anna Koeferle Dec 2015

###USAGE###
# python Plot_Log2FC_along_chr.py -c [df_Log2FC] -l [Library_File] -t

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import collections as coll
import re
import pylab
import pybedtools
import argparse
import sys

###Functions####
def load_into_df(arg):
    user_input = pd.DataFrame.from_csv(arg, header=0, sep=',', index_col=0)
    return user_input

def bedtools_intersect_gRNAs_with_library_genes(my_df, library_regions):
    my_df = my_df["log2FoldChange"]
    gRNA_list = list(my_df.index.values)
    gRNA_string = gRNA_coord_to_bed(gRNA_list)
    gRNA_bed = pybedtools.BedTool(gRNA_string, from_string=True)
    gRNA_with_gene = gRNA_bed.intersect(library_regions, f=1, wa=True, wb=True)
    gene_name = [(f[6]) for f in gRNA_with_gene]
    data={'Log2FC':my_df, 'gene':gene_name}
    df_with_target_gene=pd.DataFrame(data=data)
    return df_with_target_gene

def gRNA_coord_to_bed(gRNA_list):
    empty_list = []
    for gRNA in gRNA_list:
        lsplit = gRNA.split (':')
        chrom = lsplit[0]
        startstop = lsplit[1].split('-')
        start = startstop[0]
        stop = startstop[1]
        empty_list.append( chrom + " " + start + " " + stop + " \n")
    gRNA_string = ' '.join(empty_list)
    return gRNA_string

def find_middle_nt_of_gRNA(my_df):
    list_of_middle_pos=[]
    for i in my_df.index:
        split = i.split(':')
        start = split[1].split('-')
        middle_bp = int(start[0]) + 11
        list_of_middle_pos.append(middle_bp)
    return list_of_middle_pos

def get_list_of_genes_in_library(arg):
    EMT_lib_dict = {}
    with open(arg, 'r') as infile:
        for line in infile:
            lsplit = line.strip().split('\t')
            gRNA = lsplit[0] + ":" + lsplit[1] + "-" +lsplit[2]
            gene = lsplit[3]
            EMT_lib_dict[gRNA] = gene
    gene_names=[]
    for coord in EMT_lib_dict:
        if EMT_lib_dict[coord] not in gene_names:
            gene_names.append(EMT_lib_dict[coord])
    return gene_names

def find_distance_from_TSS(gene_name, df_this_gene, TSS_genes):
    row_TSS = TSS_genes.loc[TSS_genes['Name']==gene_name]
    if (row_TSS['strand'] == '-').bool():  #the comparison is made within the pandas dataframe, so in order to get the boolean out of the dataframe, need to apply the .bool() method
        TSS_this_gene = int(str(row_TSS['txEnd'].tolist()).strip('[]')) #select the object from dataframe and convert it to list, then convert list to string and finally convert string to int
        distance_from_tss = TSS_this_gene - df_this_gene['gRNA_central_bp']
    if (row_TSS['strand'] == '+').bool():
        TSS_this_gene = int(str(row_TSS['txStart'].tolist()).strip('[]'))
        distance_from_tss = df_this_gene['gRNA_central_bp'] - TSS_this_gene
        ###this means that no matter if the gene is on the + or - strand the distance will be negative when the gRNA is upstream and positive when downstream of the TSS
    return distance_from_tss

def plot_scatter(x_distance_from_tss, y_Log2FC, gene_name):
    fig,ax = plt.subplots() # use this to get ax object
    ax.scatter(x_distance_from_tss, y_Log2FC, c='black', edgecolor='none', s=30) # plots the scatterplot

    ax.set_xlabel('Distance from TSS (bp)', fontsize=16) # label for x axis
    ax.set_ylabel('Log2 Fold Change', fontsize=16) # label for y axis
    fig.suptitle(str(gene_name), fontsize=22) #figure title
    xdim = ax.get_xlim() # from the x object get the range of the data, returns a tuple
    xmin = xdim[0] # subset the tuple
    xmax = xdim[1]
    ydim = ax.get_ylim()
    ymin = ydim[0]
    ymax = ydim[1]
    plt.axhline(0, color='black',ls='--', lw=2)
    plt.axhline(2, color='gainsboro',ls='--', lw=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.get_yaxis().set_tick_params(direction='out', labelsize=14)
    ax.get_xaxis().set_tick_params(direction='out', labelsize=14)
    plt.locator_params(axis='x', nbins=6) #adjusts the number of ticks on the x and y axes
    plt.locator_params(axis='y', nbins=5)
    #ax.set_xlim([-1000, +500])
    plt.savefig(str(gene_name) + '_L2FC.png', bbox_inches='tight', format='png', dpi=300) ## bbox extra artist and bbox inches makes sure legend is not cut off figure
    plt.close()

##Catching user input###

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--Log2FC", help="dataframe containing Log2FC for each gRNA")
parser.add_argument("-l", "--library_file", help="tab-separated file containing library regions in bed file format")
parser.add_argument("-t", "--tss_file", help="tab-separated file containing TSS information for genes in the library")
args = parser.parse_args()

if not args.Log2FC or not args.library_file or not args.tss_file:
   sys.exit("Missing inputfile!")

print "Reading data...."
main_df = load_into_df(args.Log2FC) #loads the dataframe of DESeq2 output log2FC values
library_regions = pybedtools.BedTool(args.library_file) #loads the bedfile of regions targeted by gRNAs in the library

main_df = bedtools_intersect_gRNAs_with_library_genes(main_df, library_regions) #add a column for the target gene for each gRNA to the df
main_df['gRNA_central_bp']=find_middle_nt_of_gRNA(main_df) # add a column containing the middle chromosome coordinate of the gRNA. This is a the central nucleotide and will be used to plot gRNA position on the x axis

gene_names=get_list_of_genes_in_library(args.library_file) #get the list of genes in the library

TSS_genes=pd.DataFrame.from_csv(args.tss_file, sep=',') #load the TSS file into a dataframe

for gene_name in gene_names:
    df_this_gene = main_df.loc[main_df['gene'] == gene_name] #subset only data for single gene in list
    df_this_gene.insert(3, 'distance_from_tss', find_distance_from_TSS(gene_name, df_this_gene, TSS_genes))
    plot_scatter(df_this_gene['distance_from_tss'], df_this_gene['Log2FC'], gene_name)

print "...Done!"
