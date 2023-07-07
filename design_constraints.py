# Mia Mihailovic
# This program evaluates the impact of "pseudo" target regions (modeled as unpaired oligonucleotides) on the RBS exposure 
# (ie base pairing probability, BPP) of the long-range iRS3 system by compiling corresponding Nupack BPP predictions 
#  
# Note that Nupack was run in parallel for all sequence combinations [ (1) iRS3 only, (2) iRS3 + asRNA1, (3) iRS3 + asRNA2 and 
# (4) iRS3 + asRNA1 + asRNA2 aka "both"] on TACC using pairs: eg "pairs -T 37 -multi rna tmp/X.asRNA1"
# then standard cleanup perform was performed on all Nupack outputs (via AL's nupack_cleanup5b_args.py code) to create .clean.txt files
# output: construct sequence details with RBS base pairing probabilities (sum and st dev) for each combination for downstream user-guided 
# filtering for constructs that are fully triggered to the binding of both asRNAs, but not one or the other individually

import os
import sys
import Bio
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
import math
import statistics
import argparse
import csv

#------------------------------------------------------------------------------
# Define command line input arguments
parser = argparse.ArgumentParser(description='long range targets optimization')
parser.add_argument('-i',  
    action='store', 
    dest='i',
    required=True,
    type=str,
    help="input _iRS3.clean.txt file from NUPACK cleanup")
parser.add_argument('-asRNA1',  
    action='store', 
    dest='asRNA1',
    required=True,
    type=str,
    help="input _asRNA1.clean.txt file from NUPACK cleanup")
parser.add_argument('-asRNA2',
    action = 'store', 
    dest='asRNA2',
    required =True, 
    type=str, 
    help = "input _asRNA2.clean.txt file from NUPACK cleanup")
parser.add_argument('-both',
    action = 'store', 
    dest= 'both', 
    required =True, 
    type=str, 
    help = "input _both.clean.txt file from NUPACK cleanup")
parser.add_argument('-ind',
    action = 'store', 
    dest= 'ind', 
    required =True, 
    type=str, 
    help = "index with which to find RBS start/stop")
parser.add_argument('-o',
    action = 'store', 
    dest= 'o', 
    required =True, 
    type=str, 
    help = "output summary file name")


#------------------------------------------------------------------------------
options = parser.parse_args()

#import table containing names, RBS start, RBS end: split, append columns to list
#this csv file is required for the indices for RBS start/end for nupack result extraction
#this csv file is based off of the output from AsRNA_Probe_Design.m
names = []
all_RBS_starts = []
all_RBS_ends = []
with open('construct_details.csv','r') as csvfile:
    read_file = csv.reader(csvfile, delimiter = ",")
    for row in read_file:
        names.append(row[0])
        all_RBS_starts.append(row[7])
        all_RBS_ends.append(row[8])
csvfile.close()

#find RBS start and end for each construct, keeping in mind that indexing in python starts from 0 not 1
RBS_start = all_RBS_starts[int(options.ind)-1]
#print(RBS_start)
RBS_end = all_RBS_ends[int(options.ind)-1]
#print(RBS_end)
bppsums=[]
bppstds=[]

textfiles= [options.i,options.asRNA1,options.asRNA2,options.both]
for n in range(len(textfiles)):
    #import .clean.txt files
    linecounter = 0 
    bpp = open(textfiles[n], 'r')
    lines = bpp.readlines()
    probabilities = []
    
    #Mia edit to calculate sum of bpp
    for line in lines:
        data = line.strip().split('\t')
        probabilities.append(float(data[1]))
    bppsums.append(sum(probabilities[(int(RBS_start)-1):int(RBS_end)]))
    bppstds.append(statistics.stdev(probabilities[(int(RBS_start)-1):int(RBS_end)]))
    #print(probabilities)

print(bppsums)

construct = []
construct_check = []
iRS3_bpp_sum = []
iRS3_bpp_std = []
asRNA1_bpp_sum = []
asRNA1_bpp_std = []
asRNA2_bpp_sum = []
asRNA2_bpp_std = []
both_bpp_sum = []
both_bpp_std = []

construct.append(names[int(options.ind)-1])
construct_check.append(int(options.ind))
iRS3_bpp_sum.append(bppsums[0])
iRS3_bpp_std.append(bppstds[0])
asRNA1_bpp_sum.append(bppsums[1])
asRNA1_bpp_std.append(bppstds[1])
asRNA2_bpp_sum.append(bppsums[2])
asRNA2_bpp_std.append(bppstds[2])
both_bpp_sum.append(bppsums[3])
both_bpp_std.append(bppstds[3])

##export to a csv file here that can be concatenated with original csv
with open(options.o + "_stats"+".csv","w") as stats:
    writer = csv.writer(stats)
    writer.writerows(zip(construct,construct_check,iRS3_bpp_sum,asRNA1_bpp_sum,asRNA2_bpp_sum,both_bpp_sum,iRS3_bpp_std,asRNA1_bpp_std,asRNA2_bpp_std,both_bpp_std))
stats.close()