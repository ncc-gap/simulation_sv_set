from __future__ import print_function 
import os, sys
import subprocess
import re

def get_position(infoA, infoB):

    posA = getpos(infoA)
    posB = getpos(infoB)
    if posB == '':
        posA, posB = getposend(infoA)
    return posA, posB

def getpos(info):
    position = ''
    l_info = info.split(';')
    for v in l_info:
        if v.startswith('POS='):
            position = v.replace('POS=','')
            break
    return position

def getposend(info):
    position = ''
    endpos = ''
    l_info = info.split(';')
    for v in l_info:
        if v.startswith('POS='):
            position = v.replace('POS=','')
        if v.startswith('END='):
            endpos = v.replace('END=','')
    return position, endpos


def get_alt_seq(alt):
    ret = ""
    if alt.find("[") > -1:
        l_alt = alt.split("[")
        if len(l_alt[0]) > 1:   ret = (l_alt[0])[1:]
        elif len(l_alt[2]) > 1: ret = (l_alt[2])[:-1]
            
    elif alt.find("]") > -1:
        l_alt = alt.split("]")
        if len(l_alt[0]) > 1:   ret = (l_alt[0])[1:]
        elif len(l_alt[2]) > 1: ret = (l_alt[2])[:-1]
        
    return ret
    

def get_homlen(info):
    homlen = 0
    l_info = info.split(';')
    for v in l_info:
        if v.startswith('HOMLEN='):
            homlen = v.replace('HOMLEN=','')
    return homlen


def get_info_val(infos, key):
    ret = ""
    l_info = infos.split(';')
    for info in l_info:
        #if info.startswith(key):
        if info.startswith(key + "="):
            ret = info.split("=")[1]
            break
    return ret
    
def get_format_val(format_keys, format_vals, key):
    ret = ""
    l_format_keys = format_keys.split(':')
    l_format_vals = format_vals.split(':')
    for index, format_key in enumerate(l_format_keys):
        if format_key == key:
            ret = l_format_vals[index]
            break
    return ret


def sort_breakpoint(chr1,pos1,dir1,chr2,pos2,dir2,h_chrom_number):
    flag = sort_breakpoint_main(chr1,pos1,chr2,pos2,h_chrom_number)
    if flag:
        return chr1, pos1, dir1, chr2, pos2, dir2
    else:
        return chr2, pos2, dir2, chr1, pos1, dir1
        

def sort_breakpoint_vcfid(chr1, pos1, dir1, vcfid1, chr2, pos2, dir2, vcfid2,h_chrom_number):
    flag = sort_breakpoint_main(chr1,pos1,chr2,pos2,h_chrom_number)
    if flag:
        return chr1, pos1, dir1, vcfid1, chr2, pos2, dir2, vcfid2
    else:
        return chr2, pos2, dir2, vcfid2, chr1, pos1, dir1, vcfid1


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


def get_hash_key(F, margin):
    
    chrA, startA, endA = F[0], str(int(F[1])+int(margin)), str(int(F[2])-int(margin))
    chrB, startB, endB = F[3], str(int(F[4])+int(margin)), str(int(F[5])-int(margin))
    vcfid = F[6]
    strandA = F[8]
    strandB = F[9]
    # key 
    key = chrA+'\t'+startA+'\t'+endA+'\t'+chrB+'\t'+startB+'\t'+endB+'\t'+vcfid+'\t'+strandA+'\t'+strandB
    return key

def get_hash_sort_value(chrA, startA, strandA, chrB, startB, strandB, insert_seq, h_chrom_number):
    
    chrA, startA, strandA, chrB, startB, strandB = sort_breakpoint(chrA, startA, strandA, chrB, startB, strandB, h_chrom_number)
    val = chrA+','+startA+','+strandA+','+chrB+','+startB+','+strandB+','+insert_seq
    return val

def get_hash_key_singleend_sv(F, margin):
    
    chrA, startA, endA = F[0], str(int(F[1])+int(margin)), str(int(F[2])-int(margin))
    vcfid = F[3]
    strandA = F[5]
    # key 
    key = chrA+'\t'+startA+'\t'+endA+'\t'+vcfid+'\t'+strandA
    return key
    
def make_chrom_number_dict(vcf_header):

    #h_chrom_number = {}
    #with open(vcf_header, 'r') as hin:
    #    count = 1
    #    for line in hin:
    #        if line.startswith('##contig=<ID='): 
    #            line = line.rstrip('\n')
    #            line = line.replace("##contig=<ID=","")
    #            F = line.split(",")
    #            h_chrom_number[F[0]] = count
    #            count += 1
    #
    #return h_chrom_number

    keys = []
    with open(vcf_header, 'r') as hin:
        for line in hin:
            if line.startswith('##contig=<ID='): 
                line = line.rstrip('\n')
                line = line.replace("##contig=<ID=","")
                F = line.split(",")
                keys.append(F[0])

    h_chrom_number = {}
    count = 1
    for key in sorted(keys):
        h_chrom_number[key] = count
        count += 1
    return h_chrom_number
