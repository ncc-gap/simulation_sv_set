
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

def get_strand(sv_type):
    strand1 = []
    strand2 = []        
    if sv_type == "DEL":
        strand1.append('+')
        strand2.append('-')     
    elif sv_type == "INV":
        strand1.append('+')
        strand2.append('+')
        strand1.append('-')
        strand2.append('-')
    elif sv_type == "DUP":
        strand1.append('-')
        strand2.append('+')
    elif sv_type == "INS":
        strand1.append('+')
        strand2.append('-')     
    elif sv_type == "TRA":
        strand1.append('+')
        strand2.append('-')
        strand1.append('-')
        strand2.append('+')
    return strand1, strand2
    

infile1 = sys.argv[1]
infile2 = sys.argv[2]

max_line1 = 0
d_infile1 = {}
with open(infile1, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        d_infile1[F[4]] = line 
        max_line1 = F[4]

max_line2 = 0
d_infile2 = {}
with open(infile2, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        d_infile2[F[4]] = line 
        max_line2 = F[4]

max_line = int(max_line1) if int(max_line1) >= int(max_line2) else max_line2 

for index in range(1, max_line+1): 
    if str(index) in d_infile1 and str(index) in d_infile2:  
        line1 = d_infile1[str(index)] 
        line2 = d_infile2[str(index)] 
        F1 = line1.split('\t')
        F2 = line2.split('\t')

        chr1, start1, end1 = F1[0], F1[1], F1[2]
        chr2, start2, end2 = F2[0], F2[1], F2[2]
        sv_type = F1[3]
        strand1, strand2 = get_strand(sv_type)

        h_chrom_number = make_chrom_number_dict()
        f_sort = sort_breakpoint_main(F1[0],F1[1],F2[0],F2[1],h_chrom_number)
                    
        if not f_sort:
            chr1, start1, end1 = F2[0], F2[1], F2[2]
            chr2, start2, end2 = F1[0], F1[1], F1[2]
            
        sv_size = "---"
        if sv_type != "TRA":
            sv_size = str(int(end2)-int(end1))
        if sv_type == "INS":
            start2 = str(int(start1) + 1)
            end2 = str(int(end1) + 1)
            
        for i in range(len(strand1)):
            print(chr1+'\t'+start1+'\t'+end1+'\t'+chr2+'\t'+start2+'\t'+end2+'\t.\t.\t'+strand1[i]+'\t'+strand2[i]+'\t'+sv_type+'\t'+sv_size)

