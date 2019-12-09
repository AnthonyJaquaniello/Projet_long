import pysam as ps
import os
import sys
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from math import ceil

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

all_chr = {}

with open(path_results + "chip_seq_peaks.txt","w") as filout:
    filout.write("chrom1"+'\t'+"start1"+'\t'+"end1"+'\n')
    
for chromo in chromosoms:
    print(chromo)
    k = {}
    coverage = []
    k["chrom1"] = chromo
    for pileupcolumn in bam.pileup(chromo):
        coverage.append(pileupcolumn.n)
    peaks = find_peaks(coverage,threshold=2)
    width = peak_widths(coverage,list(peaks[0]))
    k["start1"] = list(peaks[0])
    k["end1"] = list(peaks[0] + width[0])
    all_chr[chromo] = k
    with open(path_results + "chip_seq_peaks.txt","a") as filout:
        for start,end in zip(k["start1"], k["end1"]):
            filout.write("{}\t{}\t{}\n".format(k["chrom1"],start,ceil(end)))
    print(" ok !")
