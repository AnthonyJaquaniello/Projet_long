#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 19:50:11 2019
@author: axel Kournak
In a group of genomic positions like detected loops or detected hairpin
how much are on chip seq peak from a protein?
"""
import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
from matplotlib.font_manager import FontProperties
import scipy.stats as stats

BIN = 2000
NTOT = 6036 # total number of bins for cerevisiae

name = sys.argv[3]
peaks = sys.argv[2] #file from peaks_extraction.py
#file_chip = '/home/axel/Bureau/YEAST/positions_peaks_classic_cohesin.bg2'
#bank='/home/axel/Bureau/test_chromosight/yeast/SRR8769549_mitotic'
pattern_path = sys.argv[1] #directory where file from chromosight is (without '/' at the end)

ADD_IMPR = 1

for p in ["loops", "borders", "hairpins"] :
    l1 = 0
    name_patterns = p
    # Group1
    # ChIPseq peaks 
    group1 = pd.read_table(peaks, header=None, delimiter="\t", skiprows=1) 
    group1[1] = [round(int((x+y)/2/BIN), None) for x, y in zip(group1[1], group1[2])]
    
    groupA = group1[[0, 1]]
    groupA = groupA.drop_duplicates(keep="first")
    groupA = groupA.rename(columns={0: "chrom1", 1: 'start1'})
    
    # Group2:
    # detected patterns:
    if name_patterns == "borders" :
        try:
            group2 = pd.read_table(pattern_path+'/borders_out.txt', header=0, delimiter="\t") 
            group2['start1'] = [int((x+y)/2/BIN) for x, y in zip(group2['start1'], group2['end1'])]
        except FileNotFoundError:
            pass        
    elif name_patterns == "hairpins" :
        try:
            group2 = pd.read_table(pattern_path+'/hairpins_out.txt', header=0, delimiter="\t") 
            group2['start1'] = [int((x+y)/2/BIN) for x, y in zip(group2['start1'], group2['end1'])]
        except FileNotFoundError:
            pass
    elif name_patterns == "loops" :
        try:
            group2 = pd.read_table(pattern_path+'/loops_out.txt', header=0, delimiter="\t") 
            group2['start1'] = [int((x+y)/2/BIN) for x, y in zip(group2['start1'], group2['end1'])]
            group2['start2'] = [int((x+y)/2/BIN) for x, y in zip(group2['start2'], group2['end2'])]    
            group2a = group2[['chrom1', 'start1']]
            group2b = group2[['chrom2', 'start2']]
            group2b = group2b.rename(columns={'chrom2': "chrom1", 'start2': 'start1'})
            group2 = pd.concat([group2a, group2b])
        except FileNotFoundError:
            pass
    
    for p in range(len(group2)) :  # over all patterns 
        ligne2 = group2.iloc[[p]]
        if ADD_IMPR == 1:
            ligne22 = ligne2.copy()
            ligne22['start1'] = ligne22['start1'] + 1
            ligne222 = ligne2.copy()
            ligne222['start1'] = ligne222['start1'] - 1
            frames = [ligne222, ligne2, ligne22]
            ligne2_impr = pd.concat(frames)     
            ligneB = ligne2_impr[['chrom1', 'start1']]
        else : 
            ligneB = ligne2[['chrom1', 'start1']]
            
        if len((pd.merge(groupA, ligneB)).drop_duplicates(keep="first")) >= 1 :
            l1 += 1
    
    # Overlap between both groups and statistical test:   
    l2 = len(group2)
    l3 = len(groupA)
    l4 = NTOT
    
    Mcont = [[l1,l2-l1], [l3,l4-l3]]
    oddsratio, pvalue = stats.fisher_exact(Mcont)
    print("Oddsration = {}\npvalue = {}".format(oddsratio, pvalue))
    
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
    plt.ylabel('Proportion of bins enriched in '+name)
    plt.title('Proportion of bins, pvalue~ '+str(pvalue))
    plt.xticks(ind, ('Detected '+name_patterns,'All genome'))
    plt.yticks(np.arange(0, 81, 10))
    
    plt.text(0-width/2.,cat1[0]/2.0,str(l1)+"\n"+"("+str(round(float(l1)/l2*100.,2))+"%)", fontproperties=font)
    plt.text(0-width/2.,cat1[0]+cat2[0]/2, l2-l1, fontproperties=font)
    
    plt.text(1-width/2.,cat1[1]/2.0,str(l3)+"\n"+"("+str(round(float(l3)/l4*100.,2))+"%)", fontproperties=font)
    plt.text(1-width/2.,cat1[1]+cat2[1]/2.,l4-l3, fontproperties=font)
    
    #plt.legend((pattern_path+" "+name+"_"+ ' peak', 'other'), loc=1)
    #plt.show()
    plt.savefig(pattern_path+"/proportion3_-"+name+"_"+name_patterns+'.png')
    plt.close('all')

    # pie chart:
    nb_1=  l1
    nb_2 = l2-l1

    # Plot: make a square figure and axes to plot a pieChart:
    plt.figure(1, figsize=(6,6))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    # The slices will be ordered and plotted counter-clockwise.
    labels = 'in 20% most expr bins', 'Other'
    fracs = [nb_1, nb_2]
    colors = ['gold', 'lightskyblue']
    plt.pie(fracs , labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=90);
    plt.title(' ');
    #plt.show()
    plt.savefig(pattern_path+"/"+"PieChart_events_"+name+"_"+"_"+name_patterns+'.png')
    plt.close('all')
