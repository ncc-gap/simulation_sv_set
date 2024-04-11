from __future__ import print_function
import os, sys

def filter_min_sv_size(bedpe_file, output, min_sv_size):
  
    # There should be no ambiguos range!
    # chrom1, pos1-1, pos1, chrom2, pos2-1, pos2 
    with open(output, 'w') as hout:
        with open(bedpe_file, 'r') as hin:
            for line in hin:
                
                if line.startswith("#"):
                    header = line.rstrip('\n')
                    print(header, file=hout)
                    continue
                line = line.rstrip('\n')
                F = line.split('\t')
                chr1 = F[0]
                chr2 = F[3]
                end1 = F[2]
                end2 = F[5]
                if chr1 == chr2 and  (int(end2) - int(end1)) < int(min_sv_size): continue
                print(line, file=hout)
    

def filter_scaffold(bedpe_file, output, f_grc):
    
    l_chr = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    if f_grc == True:
        l_tmp_chr = []
        for chrom in l_chr:
            l_tmp_chr.append('chr'+ chrom)
        l_chr = l_tmp_chr

    hOUT = open(output, 'w')
    with open(bedpe_file, 'r') as hin:
        for line in hin:
            
            if line.startswith("#"):
                header = line.rstrip('\n')
                print(header, file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrA = F[0]
            chrB = F[3]
            if chrA not in l_chr: continue
            if chrB not in l_chr: continue
            print(chrA+'\t'+'\t'.join(F[1:3])+'\t'+chrB+'\t'+'\t'.join(F[4:6])+'\t'+'\t'.join(F[6:]), file=hOUT)
    hOUT.close()
    

def filter_scaffold_singleend_SV(bed_file, output, f_grc):
    
    l_chr = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    if f_grc == True:
        l_tmp_chr = []
        for chrom in l_chr:
            l_tmp_chr.append('chr'+ chrom)
        l_chr = l_tmp_chr

    hOUT = open(output, 'w')
    with open(bed_file, 'r') as hin:
        for line in hin:
            
            if line.startswith("#"):
                header = line.rstrip('\n')
                print(header, file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrA = F[0]
            if chrA not in l_chr: continue
            print(chrA+'\t'+'\t'.join(F[1:3])+'\t'+'\t'.join(F[3:]), file=hOUT)
    hOUT.close()

def filter_bedpe_main(args):

    filter_scaffold(args.in_bedpe, args.output, args.grc)
