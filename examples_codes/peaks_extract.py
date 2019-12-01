import pysam as ps
import os
import sys
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

name_correcting = {'chr1':'NC_001133.9','chr2':'NC_001134.8','chr3':'NC_001135.5',
                   'chr4':'NC_001136.10','chr5':'NC_001137.3','chr6':'NC_001138.5',
                   'chr7':'NC_001139.9','chr8':'NC_001140.6','chr9':'NC_001141.2',
                   'chr10':'NC_001142.9','chr11':'NC_001143.9','chr12':'NC_001144.5',
                   'chr13':'NC_001145.3','chr14':'NC_001146.8','chr15':'NC_001147.6',
                   'chr16':'NC_001148.4'}

path_results = sys.argv[2]

name_conversion = {}
for key in name_correcting:
    name_conversion[name_correcting[key]] = key
    
#sys.argv[1] = path to bamfile
bam = ps.AlignmentFile(sys.argv[1],'rb')
chromosoms = list(name_correcting.values())
all_chr = {} #dict to dict : key = one chr, value = dict (see below)

for chromo in chromosoms:
    coverage = {} #dict: key = pos, value = coverage
    for pileupcolumn in bam.pileup(chromo):
        coverage[pileupcolumn.pos] = pileupcolumn.n
    all_chr[chromo] = coverage
    vect = list(all_chr[chromo].values()) #coverage values
    peaks = find_peaks(vect,width= 1,threshold=3)   
    with open(path_results+name_conversion[chromo]+'_peaks.txt','w') as filout:
        for peak in list(peaks[0]):
            filout.write(str(peak)+"\n")
