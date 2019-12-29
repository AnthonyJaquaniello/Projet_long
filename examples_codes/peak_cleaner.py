import sys
import pandas as pd

peak_real = sys.argv[1] #ChIP-Seq peaks with IP
peak_input = sys.argv[2] #ChIP-Seq peaks without IP
path = sys.argv[3]

chrom_real = []
start_real = []
end_real = []
chrom_input = []
start_input = []
end_input = []

with open(peak_real, 'r') as real:
    real.readline() #pas de headers
    for line in real:
        chrom_real.append(line.split()[0])
        start_real.append(int(line.split()[1]))
        end_real.append(int(line.split()[2]))
        
with open(peak_input, 'r') as inputt:
    inputt.readline()
    for line in inputt:
        chrom_input.append(line.split()[0])
        start_input.append(int(line.split()[1]))
        end_input.append(int(line.split()[2]))

pos_real = [int((x+y)/2) for x, y in zip(start_real, end_real)]
pos_input = [int((x+y)/2) for x, y in zip(start_input, end_input)]

d_real = {"Name":chrom_real,"Moy":pos_real,
	  "Start":start_real, "End":end_real}
d_input = {"Name":chrom_input,"Moy":pos_input,
	   "Start":start_input, "End":end_input}

real = pd.DataFrame(d_real)
inputt = pd.DataFrame(d_input)

left_real = real.merge(right=inputt, left_on=['Name','Moy'],
		       right_on=['Name','Moy'],how='left')

true_real = left_real[left_real["Start_y"].isna()]
true_real = true_real.drop("Start_y",axis=1)
true_real = true_real.drop("End_y", axis=1)
true_real = true_real.drop("Moy", axis=1)
true_real.to_csv(path_or_buf=path,sep="\t",
		 header=["chrom1","start1","end1"], index= False)


