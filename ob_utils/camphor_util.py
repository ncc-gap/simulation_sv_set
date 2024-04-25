import os, sys
import subprocess, shutil
import utils
from svtools.vcftobedpe import run_vcf2bedpe 
from filter_bedpe import filter_scaffold
import pysam
import re

def __get_INFO(F):
    info = {}
    for cell in F[18].split(";"):
        if "=" in cell:
            (key, value) = cell.split("=")
            info[key] = value
        else:
            info[cell] = True
    return info

def camphorSVtoBedpe(input_vcf, output, patched_vcf, f_grc, filter_scaffold_option, bcf_filter_option, debug):

    out_pref, ext = os.path.splitext(output)

    hOUT = open(patched_vcf, 'w')
    svtype_list = []
    with open(input_vcf, 'r') as hin:
        print("##fileformat=VCFv4.2", file=hOUT)
        print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">', file=hOUT)
        print('##INFO=<ID=LEN,Number=1,Type=Integer,Description="">', file=hOUT)
        print('##contig=<ID=chr1,length=248956422>', file=hOUT)
        print('##contig=<ID=chr2,length=242193529>', file=hOUT)
        print('##contig=<ID=chr3,length=198295559>', file=hOUT)
        print('##contig=<ID=chr4,length=190214555>', file=hOUT)
        print('##contig=<ID=chr5,length=181538259>', file=hOUT)
        print('##contig=<ID=chr6,length=170805979>', file=hOUT)
        print('##contig=<ID=chr7,length=159345973>', file=hOUT)
        print('##contig=<ID=chr8,length=145138636>', file=hOUT)
        print('##contig=<ID=chr9,length=138394717>', file=hOUT)
        print('##contig=<ID=chr10,length=133797422>', file=hOUT)
        print('##contig=<ID=chr11,length=135086622>', file=hOUT)
        print('##contig=<ID=chr12,length=133275309>', file=hOUT)
        print('##contig=<ID=chr13,length=114364328>', file=hOUT)
        print('##contig=<ID=chr14,length=107043718>', file=hOUT)
        print('##contig=<ID=chr15,length=101991189>', file=hOUT)
        print('##contig=<ID=chr16,length=90338345>', file=hOUT)
        print('##contig=<ID=chr17,length=83257441>', file=hOUT)
        print('##contig=<ID=chr18,length=80373285>', file=hOUT)
        print('##contig=<ID=chr19,length=58617616>', file=hOUT)
        print('##contig=<ID=chr20,length=64444167>', file=hOUT)
        print('##contig=<ID=chr21,length=46709983>', file=hOUT)
        print('##contig=<ID=chr22,length=50818468>', file=hOUT)
        print('##contig=<ID=chrX,length=156040895>', file=hOUT)
        print('##contig=<ID=chrY,length=57227415>', file=hOUT)

        for line in hin:
            line = line.rstrip('\n')
            if line.startswith("#"):
                print(line, file=hOUT)
                continue
            F = line.split('\t')
            svtype = F[4].replace("<", "").replace(">", "")
            print(line + ";SVTYPE=" + svtype, file=hOUT)
            if not svtype in svtype_list:
                svtype_list.append(svtype)
        print("SVTYPE is " + ",".join(sorted(svtype_list)))
    hOUT.close()

    if bcf_filter_option != "":
        subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", "-f", bcf_filter_option, patched_vcf])
    else:
        subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", patched_vcf])
        
    run_vcf2bedpe(out_pref + ".tmp1.vcf", out_pref + ".tmp1.bedpe")

    with open(out_pref + ".tmp1.bedpe", 'r') as hin:
        with open(out_pref + ".tmp2.bedpe", 'w') as hout:
            for line in hin:
                line = line.rstrip('\n')
                if line.startswith("#"):
                    print(line, file=hout)
                    continue
                F = line.split("\t")
                info = __get_INFO(F)
                if "CHR2" in info:
                    F[3] = info["CHR2"]
                F[8] = "*"
                F[9] = "*"
                print('\t'.join(F), file=hout)

    if filter_scaffold_option:
        filter_scaffold(out_pref + ".tmp2.bedpe", out_pref + ".tmp3.bedpe", f_grc)
    else:
        shutil.copyfile(out_pref + '.tmp2.bedpe', out_pref + '.tmp3.bedpe')

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".tmp3.bedpe"], stdout = hOUT)
    hOUT.close()

    if not debug:
        os.remove(out_pref + ".tmp1.vcf")
        os.remove(out_pref + ".tmp1.bedpe")
        os.remove(out_pref + ".tmp2.bedpe")
        os.remove(out_pref + ".tmp3.bedpe")

def filt_clustered_rearrangement2(input_file, output_file, min_tumor_support_read, min_sv_length, h_chrom_number):

    hout = open(output_file, 'w')
    control_support_read = "---"
    tumor_ref_read = "---"
    insert_seq = "---"
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            
            tchr1, tstart1, tend1, tchr2, tstart2, tend2 = F[0], str(int(F[1])-1), F[2], F[3], str(int(F[4])-1), F[5]
            tdir1, tdir2 = F[8], F[9]
            sv_type = F[10]
            
            alt = F[14]
            info1, info2 = F[18], F[19]
            
            format_keys = ""
            format_vals = ""
            tumor_support_read = utils.get_info_val(info1, "RE")
            if int(tumor_support_read) < min_tumor_support_read: continue
            control_ref_read = utils.get_info_val(info1, "ND")
            
            sort_flag = utils.sort_breakpoint_main(tchr1,tstart1,tchr2,tstart2,h_chrom_number)
            if not sort_flag:
                tchr1, tstart1, tend1, tdir1, info1, tchr2, tstart2, tend2, tdir2, info2 = tchr2, tstart2, tend2, tdir2, info2, tchr1, tstart1, tend1, tdir1, info1

            print("\t".join([tchr1, tend1, tdir1, tchr2, tend2, tdir2, insert_seq, tumor_ref_read, tumor_support_read, control_ref_read, control_support_read, sv_type]), file = hout)

    hout.close()

def camphorSVtoBedpe_main(args):
    
    in_tumor_sv = args.in_tumor_sv
    f_grc = args.f_grc
    filter_scaffold_option = args.filter_scaffold_option
    bcf_filter_option = args.bcf_filter_option
    min_tumor_support_read = args.min_tumor_support_read
    min_sv_length = args.min_sv_length
    output = args.output
    debug = args.debug

    output_prefix, ext = os.path.splitext(output)
    
    camphorSVtoBedpe(in_tumor_sv, output_prefix+'.camphor_tumor_PASS.bedpe', output_prefix+'.camphor_patched.vcf', f_grc, filter_scaffold_option, bcf_filter_option, debug)

    bcftools_command = ["bcftools", "view", "-h", output_prefix+'.camphor_patched.vcf', "-o", output_prefix +'.camphor_patched.vcf.header']
    subprocess.check_call(bcftools_command)
    h_chrom_number = utils.make_chrom_number_dict(output_prefix +'.camphor_patched.vcf.header')
    filt_clustered_rearrangement2(output_prefix+'.camphor_tumor_PASS.bedpe', output_prefix+'.camphor_filtered.txt',
    min_tumor_support_read, min_sv_length, h_chrom_number)

    with open(output_prefix + ".camphor_sorted.txt", 'w') as hout:
        subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", "-k3,3", "-k6,6", output_prefix + ".camphor_filtered.txt"],  stdout = hout)

    with open(output, 'w') as hout:
        l_header = ["Chr_1","Pos_1","Dir_1","Chr_2","Pos_2","Dir_2","Inserted_Seq","Checked_Read_Num_Tumor","Supporting_Read_Num_Tumor","Checked_Read_Num_Control","Supporting_Read_Num_Control","Sv_Type"]
        print("\t".join(l_header), file=hout)
        with open(output_prefix + ".camphor_sorted.txt", 'r') as hin:
            for line in hin:
                print(line.rstrip('\n'), file=hout)
            
    if not debug:
        os.remove(output_prefix +'.camphor_tumor_PASS.bedpe')
        os.remove(output_prefix +'.camphor_patched.vcf')
        os.remove(output_prefix +'.camphor_patched.vcf.header')
        os.remove(output_prefix +'.camphor_filtered.txt')
        os.remove(output_prefix +'.camphor_sorted.txt')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog = "ob_utils")
    parser.add_argument("--in_tumor_sv", help = "the vcf format file", type = str, required=True)
    parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
    parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
    parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
    parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
    parser.add_argument("--min_tumor_support_read", help = "minimum tumor support reads", type = int, default = 3)
    parser.add_argument("--min_sv_length", help = "minimum sv length", type = int, default = 1)
    parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")
    args = parser.parse_args()

    camphorSVtoBedpe_main(args)
