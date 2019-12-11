#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 19:50:11 2019
@author: axel Kournak
In a group of genomic positions like detected loops or detected hairpin
how much are on chip seq peak from a protein?
"""
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.font_manager import FontProperties
import scipy.stats as stats

BIN = 2000 #10 000 for human
NTOT = 6037 #283188 total number of bins for human
pattern = sys.argv[1] #from chromosight
peaks = sys.argv[2] #from peaks_extraction.py
name = sys.argv[3]

# Group1
# ChIPseq peaks 
group1 = pd.read_table(peaks, header=None, delimiter="\t", skiprows=1)
group1[1] = [(int(x)+int(y))/2/BIN for x, y in zip(group1[1], group1[2])]

groupA = group1[[0,1]]
groupA = groupA.drop_duplicates(keep="first")
groupA = groupA.rename(columns={0: "chrom1", 1: 'start1'})

# Group2:
# detected patterns
group2 = pd.read_table(pattern, header=0, delimiter="\t") 
group2['start1'] = [int((x+y)/2/BIN) for x, y in zip(group2['start1'], group2['end1'])]
#mean between start position and end position(normalised by BIN)

# Adding some imprecision: 
group22 = group2.copy()
group22['start1'] = group2['start1']+1
group222 = group2.copy()
group222['start1'] = group2['start1']-1

frames = [group222, group2, group22]
group2_impr = pd.concat(frames)

groupB = group2[['chrom1', 'start1']]
groupB = groupB.drop_duplicates(keep="first")

# Overlap between both groups and statistical test:

l1=len((pd.merge(groupA, groupB)).drop_duplicates(keep="first"))
l2=len(groupB)
l3=len(groupA)
l4=NTOT

Mcont = [[l1, l2-l1], [l3-l1, l4-l3]]
oddsratio, pvalue = stats.fisher_exact(Mcont)

# Plots: 
font0 = FontProperties()
font = font0.copy()
font.set_weight('bold')

N = 2
cat1 = (float(l1)/l2*100., float(l3)/l4*100.)
cat2 = (l2/l2*100.-float(l1)/l2*100., l4/l4*100.-float(l3)/l4*100.)
ind = np.arange(N)    # the x locations for the groups
width = 0.3       # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind, cat1, width)
p2 = plt.bar(ind, cat2, width, bottom=cat1)

plt.xlim(-1,2.5)
plt.ylabel('Proportion of bins enriched in protein ')
plt.title('Proportion of bins enriched in protein, pvalue~ '+str(pvalue))
plt.xticks(ind, ('Detected hairpins', 'All genome'))
plt.yticks(np.arange(0, 81, 10))

plt.text(0-width/2.,cat1[0]/2.0,str(l1)+"\n"+"("+str(round( float(l1)/l2*100.,2))+"%)", fontproperties=font)
plt.text(0-width/2.,cat1[0]+cat2[0]/2, l2-l1, fontproperties=font)

plt.text(1-width/2.,cat1[1]/2.0,str(l3)+"\n"+"("+str(round( float(l3)/l4*100.,2))+"%)", fontproperties=font)
plt.text(1-width/2.,cat1[1]+cat2[1]/2.,l4-l3, fontproperties=font)

plt.legend(('protein peak', 'other'), loc=1)
plt.savefig('protein_enrichment_{}.png'.format(name))
plt.close('all')


