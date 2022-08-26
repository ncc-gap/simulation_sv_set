#! /usr/bin/env python
import sys, gzip, math
import os, glob, re

input_dir = sys.argv[1]
output = sys.argv[2]

l_infile = glob.glob(input_dir+'_support_*/*goldendata.txt')
with open(output, "w") as hout:
    
    print("\tTP\tFP\tFN\tSupport",file=hout)

    for input_file in l_infile:
        input_file_prefix, ext = os.path.splitext(input_file)
        dirname = os.path.dirname(input_file)
        tmp_name, support_read = dirname.split("_support_")

        l_true_positive = []
        l_false_positive = []
        l_false_negative = []

        with open(input_file, 'r') as hin:
            for line in hin:
                if line.startswith("Chr_1"):continue
                line = line.rstrip('\n')
                F = line.split('\t')
                key = ",".join(F[0]+'\t'+F[1]+'\t'+F[2]+'\t'+F[3]+'\t'+F[4]+'\t'+F[5])
                if F[12] != "":
                    l_true_positive.append(key)
                else:
                    l_false_positive.append(key)
                
        with open(input_file_prefix + ".false_negative.txt", "r") as hin:
            for line in hin:
                line = line.rstrip('\n')
                F = line.split('\t')
                key = ",".join(F[0]+'\t'+F[2]+'\t'+F[3]+'\t'+F[5]+'\t'+F[8]+'\t'+F[9])
                l_false_negative.append(key)

        result_name = re.findall(r'_(DP\d+)_(TP\d+)_', input_file)
        DP, TP = result_name[0]
        tp = str(len(set(l_true_positive)))
        fp = str(len(set(l_false_positive)))
        fn = str(len(set(l_false_negative)))
        print(DP+"_"+TP+'\t'+tp+'\t'+fp+'\t'+fn + '\t'+ support_read, file=hout)

'''
d_true_positive_count = {}
d_false_positive_count = {}
d_false_negative_count = {}

for key in d_true_positive:
    sv_type = d_true_positive[key]
    if sv_type not in d_true_positive_count: 
        d_true_positive_count[sv_type] = 0
    d_true_positive_count[sv_type] += 1        

for key in d_false_positive:
    sv_type = d_false_positive[key]
    if sv_type not in d_false_positive_count: 
        d_false_positive_count[sv_type] = 0
    d_false_positive_count[sv_type] += 1        

for key in d_false_negative:
    sv_type = d_false_negative[key]
    if sv_type not in d_false_negative_count: 
        d_false_negative_count[sv_type] = 0
    d_false_negative_count[sv_type] += 1 
    
print(d_true_positive_count)
print(d_false_positive_count)
print(d_false_negative_count)
'''


