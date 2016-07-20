#!/usr/bin/python
import os
import subprocess
import csv

###
# Compares IMSEQ output of pooled fastq file from ImmunoSeq enrichment protocol to
# processed output results tables provided from Adaptive
###

input_imseq = "/scratch/adaptive/adaptive_out.tsv"
input_adaptive = "/scratch/adaptive/FromImmunoSeq/"
output = "/scratch/adaptive/output/"

files = os.listdir(input_adaptive)

with open(input_adaptive + files[0], "rb") as f:
	reader = csv.reader(f, delimiter='\t')
	adaptive = list(reader)

X4015B = adaptive
X4015B_aa = [item[1] for item in adaptive]

with open(input_adaptive + files[1], "rb") as f:
	reader = csv.reader(f, delimiter='\t')
	adaptive = list(reader)

X1285F2 = adaptive
X1285F2_aa = [item[1] for item in adaptive]

with open(input_adaptive + files[2], "rb") as f:
	reader = csv.reader(f, delimiter='\t')
	adaptive = list(reader)

X1285B = adaptive
X1285B_aa = [item[1] for item in adaptive]

with open(input_adaptive + files[3], "rb") as f:
	reader = csv.reader(f, delimiter='\t')
	adaptive = list(reader)

X1304F = adaptive
X1304F_aa = [item[1] for item in adaptive]

with open(input_adaptive + files[4], "rb") as f:
	reader = csv.reader(f, delimiter='\t')
	adaptive = list(reader)

X1285F1 = adaptive
X1285F1_aa = [item[1] for item in adaptive]

with open(input_adaptive + files[5], "rb") as f:
	reader = csv.reader(f, delimiter='\t')
	adaptive = list(reader)

X1304B = adaptive
X1304B_aa = [item[1] for item in adaptive]

with open(input_adaptive + files[6], "rb") as f:
	reader = csv.reader(f, delimiter='\t')
	adaptive = list(reader)

X4015F = adaptive
X4015F_aa = [item[1] for item in adaptive]

with open(input_imseq, "rb") as g:
	reader = csv.reader(g, delimiter='\t')
	imseq = list(reader)

imseq_aa = [item[10] for item in imseq]
imseq_aa = imseq_aa[1:]

counts = dict()
for i in imseq_aa:
  counts[i] = counts.get(i, 0) + 1

freqs = list()
# aa, imseq count, freq for each Adaptive
freqs.append(["aa", "adaptive", "4015B", "1285F2", "1285B", "1304F", "1285F1", "1304B", "4015F"])
for i in counts:
	freq_row = ["",0,0,0,0,0,0,0,0]
	freq_row[0] = i
	freq_row[1] = float(counts[i])/float(len(imseq_aa))
	if (i in X4015B_aa):
		freq_row[2] = float(X4015B[X4015B_aa.index(i)][3])
	if (i in X1285F2_aa):
		freq_row[3] = float(X1285F2[X1285F2_aa.index(i)][3])
	if (i in X1285B_aa):
		freq_row[4] = float(X1285B[X1285B_aa.index(i)][3])		
	if (i in X1304F_aa):
		freq_row[5] = float(X1304F[X1304F_aa.index(i)][3])		
	if (i in X1285F1_aa):
		freq_row[6] = float(X1285F1[X1285F1_aa.index(i)][3])		
	if (i in X1304B_aa):
		freq_row[7] = float(X1304B[X1304B_aa.index(i)][3])
	if (i in X4015F_aa):
		freq_row[8] =float( X4015F[X4015F_aa.index(i)][3])		
	freqs.append(freq_row)

with open("adaptive_imseq_comparison.tsv", "wb") as g:
	writer = csv.writer(g)
	writer.writerows(freqs)	
 	