

import sys

in_file = sys.argv[1] # GRCh38.d1.vd1.sizes
out_pref = sys.argv[2]

h_chrom_size = {}
all_chrom_size = 0
with open(in_file,'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        F = line.split('\t')
        
        h_chrom_size[F[0]] = int(F[1])
        all_chrom_size += int(F[1])
        if F[0] == 'chrY':
            break
    
DUPLICATION_number = 2500
INDEL_number =  5000
INVERSION_number = 100
INV_del_number = 50
INV_dup_number = 50

l_chrom = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
for i in l_chrom:
    percentage = h_chrom_size["chr"+i] / all_chrom_size 
    dup_mumber = round(DUPLICATION_number * percentage)
    indel_mumber = round(INDEL_number * percentage)
    inv_mumber = round(INVERSION_number * percentage)
    invdel_mumber = round(INV_del_number * percentage)
    invdup_mumber = round(INV_dup_number * percentage)
    
    with open(out_pref+"_chr"+i,'w') as hout:
        dup_str = "DUPLICATION_minimum_length: 100\nDUPLICATION_maximum_length: 10000\nDUPLICATION_number: %d\n" % dup_mumber
        indel_str = "INDEL_minimum_length: 100\nINDEL_maximum_length: 10000\nINDEL_number: %d\n" % indel_mumber
        trans_str = "TRANSLOCATION_minimum_length: 0\nTRANSLOCATION_maximum_length: 0\nTRANSLOCATION_number: 0\n"
        inv_str = "INVERSION_minimum_length: 10\nINVERSION_maximum_length: 10000\nINVERSION_number: %d\n" % inv_mumber
        invdel_str = "INV_del_minimum_length: 100\nINV_del_maximum_length: 10000\nINV_del_number: %d\n" % invdel_mumber
        invdup_str = "INV_dup_minimum_length: 100\nINV_dup_maximum_length: 10000\nINV_dup_number: %d" % invdup_mumber
        print(dup_str+indel_str+trans_str+inv_str+invdel_str+invdup_str, file=hout)

