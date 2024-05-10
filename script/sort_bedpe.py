#! /usr/bin/env python
import sys

def make_chrom_number_dict():

    l_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    h_chrom_number = {}
    count = 1
    for chrom in l_chr:
        h_chrom_number["chr"+chrom] = count
        count += 1
    return h_chrom_number

def sort_breakpoint_main(chr1,pos1,chr2,pos2,h_chrom_number):
    
    # The result of SvABA has error records.:
    if "" in [chr1, pos1, chr2, pos2]:
        return True
    
    if chr1 == chr2 and int(pos1) > int(pos2):
        return False
    elif int(h_chrom_number[chr1]) > int(h_chrom_number[chr2]):
        return False
    else:
        return True
        

input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] == "Chr_1":
            print('\t'.join(F))
            continue
        tchr1, tpos1, tdir1, tchr2, tpos2, tdir2 = F[0], F[1], F[2], F[3], F[4], F[5]
        if tchr1.startswith("chrUn"): continue
        if tchr2.startswith("chrUn"): continue
        if tchr1.startswith("chrM"): continue
        if tchr2.startswith("chrM"): continue
        if tchr1.endswith("random"): continue
        if tchr2.endswith("random"): continue

        h_chrom_number = make_chrom_number_dict()
        f_sort = sort_breakpoint_main(tchr1, tpos1,tchr2,tpos2,h_chrom_number)
        if not f_sort:
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2 = tchr2, tpos2, tdir2,tchr1, tpos1, tdir1
        print('\t'.join([tchr1, tpos1, tdir1, tchr2, tpos2, tdir2])+'\t'+'\t'.join(F[6:]))
        
