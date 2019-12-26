# -*- coding: utf-8 -*-
"""
@author: axel KournaK 
To vizualise detected loops! 
"""
import time
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import random as rd
import sys
import os
import scn #matrix normalisation (remove bins for the normalization)
import ice_mirny3 #matrix normalisation
import scipy
import scipy.ndimage
import scipy.io as sio
import distance_law_human #to compute a distance law plot from a matrix
import hicstuff as hcs #toolbox for some HIC data processing
import numpy as np
import json

bank = sys.argv[1] #path toward .npz objects
#we have one .npz object for each couple of chr
pattern = sys.argv[2] #file of pattern
location = sys.argv[3]
filename = "detected_loops"

# Input of loops detected by chromosight: 
df = pd.read_table(pattern, header=None, delimiter="\t", skiprows=1)
bin_matrice = 2000 # size of bin of matrice (in bp)

# Loading and normalisation of contact matrices:
list_all_chrms= np.unique(df[0])
print(list_all_chrms)

matscn = {} #dict of matrix (value), key = chr name
th_sum = {}
indices_matrice = {}

name_correcting = {'chr1':'NC_001133.9','chr2':'NC_001134.8','chr3':'NC_001135.5',
                   'chr4':'NC_001136.10','chr5':'NC_001137.3','chr6':'NC_001138.5',
                   'chr7':'NC_001139.9','chr8':'NC_001140.6','chr9':'NC_001141.2',
                   'chr10':'NC_001142.9','chr11':'NC_001143.9','chr12':'NC_001144.5',
                   'chr13':'NC_001145.3','chr14':'NC_001146.8','chr15':'NC_001147.6',
                   'chr16':'NC_001148.4'}

for c in list_all_chrms : 
    print(c)
    #[below]Load a sparse matrix from a file using .npz format
    try:
        matraw = scipy.sparse.load_npz(bank+"/"+c+"_"+c+"_sparse_matrice_"+str(bin_matrice)+".txt.npz") #matrice brute
    except:
        print("Please check chromosomes names corresponding between pattern file and .npz file, continue...")
        matraw = scipy.sparse.load_npz(bank+"/"+name_correcting[c]+"_"+name_correcting[c]+"_sparse_matrice_"+str(bin_matrice)+".txt.npz")
    matraw = matraw + np.transpose(matraw) #a quoi ca sert d'additionner une matrice et sa transposée ?
    #C'est par rapport à la symétrie des cartes de contact... S'il y a contact entre le chr1 et le chr5 alors on veut additionner
    #le contact entre (chr1,chr5) et (chr5,chr1) car c'est la même chose ! ! !
    matscn[c] = hcs.normalize_sparse(matraw) 
    matscn[c] = matscn[c].tolil() #converting in "liste chaînée"
    try:
        matraw = scipy.sparse.load_npz(bank+"/"+c+"_"+c+"_sparse_matrice_"+str(bin_matrice)+".txt.npz")
    except:
        matraw = scipy.sparse.load_npz(bank+"/"+name_correcting[c]+"_"+name_correcting[c]+"_sparse_matrice_"+str(bin_matrice)+".txt.npz")
    matraw = matraw + np.transpose(matraw)
    th_sum[c] = np.median(np.array(matraw.sum(axis=0))) - 1.* np.std(matraw.sum(axis=0))
    indices_matrice[c] = np.where(matraw.sum(axis=0) > th_sum[c] )

# Looking for the patterns:
n_pos_set = 0
cj=0

for chr1 in list_all_chrms :
    cj+=1
    n1 = matscn[chr1].shape[0]
    print("Number of bins and after filtering poor interacting bins:")
    print(n1)
     
    mat = matscn[chr1].toarray()
    plt.figure(cj)
    plt.imshow(mat**0.25, interpolation="none", cmap="afmhot_r")
    plt.title(chr1+"\n"+bank)
    plt.xlabel("Position along the chromosome , bin "+str(bin_matrice))
    plt.savefig("{}/map_{}.png".format(location,chr1))
         
    pos_set = df.loc[(df[0] == chr1)]
    pos_set = np.array(pos_set)
    n_pos_set = n_pos_set + pos_set.shape[0]
    ns =0
    #print(pos_set)
    #print("pos_set.shape:{}".format(pos_set.shape))
    #for i in range(pos_set.shape[0]):
        #try:
            #site1 = int(pos_set[i,1])
            #site2 = int(pos_set[i,4])
        #except IndexError:
            #site1 = int(pos_set[i,1])
            #site2 = int(pos_set[i,4])
        #site1 = int(site1 / bin_matrice)
        #site2 = int(site2 / bin_matrice)
        #plot(site1, site2, "", color="blue")
        #plt.scatter(site1, site2, s=20, facecolors='none', edgecolors='b')
        #ns +=1
        #pi =0

#pdf = matplotlib.backends.backend_pdf.PdfPages(filename+"_"+str( int(area*bin_matrice/1000) )+".pdf")
#for fig in range(1, plt.figure().number): ## will open an empty extra figure :(
    #pdf.savefig(fig)
#pdf.close()
     
#plt.close("all")

print("ALL the experiements are finished!")

