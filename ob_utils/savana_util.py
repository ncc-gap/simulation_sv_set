import os, sys
import subprocess, shutil
import utils
from svtools.vcftobedpe import run_vcf2bedpe 
from filter_bedpe import filter_scaffold
import pysam
import re

def savanaSVtoBedpe(input_vcf, output, f_grc, filter_scaffold_option, bcf_filter_option):

    out_pref, ext = os.path.splitext(output)

    if bcf_filter_option != "":
        subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", "-f", bcf_filter_option, input_vcf])
    else:
        subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", input_vcf])
    
    # add END tag
    with open(out_pref + ".tmp1.vcf") as hin, open(out_pref + ".tmp2.vcf", "w") as hout:
        for row in hin:
            if row.startswith("#"):
                hout.write(row)
                if row.startswith("##INFO=<ID=END_EVENT_SIZE_MEAN"):
                    hout.write('##INFO=<ID=END,Number=1,Type=Float,Description="End">\n')
                continue
            F = row.split("\t")
            pos = F[1]
            alt = F[4]
            info = F[7]
            end = ""
            if "SVTYPE=BND" in info:
                end = alt.split(":")[1].split("[")[0].split("]")[0]
            elif "SVTYPE=INS" in info:
                end = pos
            else:
                raise Exception("Undefine %s" % (info.split(";")[0]))
            F[7] += ";END=" + end
            hout.write("\t".join(F))

    run_vcf2bedpe(out_pref + ".tmp2.vcf", out_pref + ".tmp2.bedpe")
    
    if filter_scaffold_option:
        filter_scaffold(out_pref + ".tmp2.bedpe", out_pref + ".tmp3.bedpe", f_grc)
    else:
        shutil.copyfile(out_pref + '.tmp2.bedpe', out_pref + '.tmp3.bedpe')

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".tmp3.bedpe"], stdout = hOUT)
    hOUT.close()

    os.remove(out_pref + ".tmp1.vcf")
    os.remove(out_pref + ".tmp2.vcf")
    os.remove(out_pref + ".tmp2.bedpe")
    os.remove(out_pref + ".tmp3.bedpe")

def repair_dup_strand(bedpe_file, output):
    
    hOUT = open(output, 'w')
    with open(bedpe_file, 'r') as hin:
        for line in hin:
            
            if line.startswith("#"):
                header = line.rstrip('\n')
                print(header, file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')

            sv_type = F[10]
            alt = F[14]

            if sv_type == "DEL":
                F[8] = '+'
                F[9] = '-'
            elif sv_type == "INV":
                F[8] = '*'
                F[9] = '*'
            elif sv_type == "DUP":
                F[8] = '-'
                F[9] = '+'
            elif sv_type == "INS":
                F[8] = '+'
                F[9] = '-'

                alt = ""
                for item in F[18].split(";"):
                    if item.startswith("SVLEN="):
                        sv_len = int(item.split("=")[1])
                        alt = "N" * sv_len
                if alt != "":
                    F[14] = alt

            #elif sv_type == "BND":
            #    if alt.startswith('N'):
            #        F[8] = '+'
            #    else:
            #        F[8] = '-'
            #    alt_seq = re.findall(r'([][])(.+?)([][])', alt)
            #    if alt_seq[0][0].startswith(']'):
            #        F[9] = '+'
            #    else:
            #        F[9] = '-'
            
            print('\t'.join(F), file=hOUT)

    hOUT.close()  
    

def filt_clustered_rearrangement2(input_file, output_file, min_tumor_support_read, min_sv_length, h_chrom_number):

    hout = open(output_file, 'w')
    tumor_ref_read = "---"
    control_ref_read = "---"
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr1, tstart1, tend1, tchr2, tstart2, tend2 = F[0], str(int(F[1])-1), F[2], F[3], str(int(F[4])-1), F[5]
            tdir1, tdir2 = F[8], F[9]
            sv_type = F[10]
            alt = F[14]
            info1, info2 = F[18], F[19]
            format_keys, format_vals = F[20], F[21]
            tumor_support_read = utils.get_info_val(info1, "TUMOUR_SUPPORT")
            control_support_read = utils.get_info_val(info1, "NORMAL_SUPPORT")
            insert_seq = alt if sv_type == "INS" else "---"
            if int(tumor_support_read) < min_tumor_support_read: continue

            sort_flag = utils.sort_breakpoint_main(tchr1,tstart1,tchr2,tstart2,h_chrom_number)
            if not sort_flag:
                tchr1, tstart1, tend1, tdir1, info1, tchr2, tstart2, tend2, tdir2, info2 = tchr2, tstart2, tend2, tdir2, info2, tchr1, tstart1, tend1, tdir1, info1

            print("\t".join([tchr1, tend1, tdir1, tchr2, tend2, tdir2, insert_seq, tumor_ref_read, tumor_support_read, str(control_ref_read), str(control_support_read),sv_type]), file = hout)

    hout.close()

def savanaSVtoBedpe_main(args):
    
    in_tumor_sv = args.in_tumor_sv
    f_grc = args.f_grc
    filter_scaffold_option = args.filter_scaffold_option
    bcf_filter_option = args.bcf_filter_option
    min_tumor_support_read = args.min_tumor_support_read
    min_sv_length = args.min_sv_length
    output = args.output
    debug = args.debug

    output_prefix, ext = os.path.splitext(output)
    
    savanaSVtoBedpe(in_tumor_sv, output_prefix+'.savana_tumor_PASS.bedpe', f_grc, filter_scaffold_option, bcf_filter_option)

    repair_dup_strand(output_prefix+'.savana_tumor_PASS.bedpe', output_prefix+'.savana_tumor_repaired.bedpe')

    bcftools_command = ["bcftools", "view", "-h", in_tumor_sv, "-o", output_prefix +'.savana.vcf.header']
    subprocess.check_call(bcftools_command)
    h_chrom_number = utils.make_chrom_number_dict(output_prefix +'.savana.vcf.header')
    filt_clustered_rearrangement2(output_prefix+'.savana_tumor_repaired.bedpe', output_prefix+'.savana_filtered.txt', 
    min_tumor_support_read, min_sv_length, h_chrom_number)

    with open(output_prefix + ".savana_sorted.txt", 'w') as hout:
        subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", "-k3,3", "-k6,6", output_prefix + ".savana_filtered.txt"],  stdout = hout)

    with open(output, 'w') as hout:
        l_header = ["Chr_1","Pos_1","Dir_1","Chr_2","Pos_2","Dir_2","Inserted_Seq","Checked_Read_Num_Tumor","Supporting_Read_Num_Tumor","Checked_Read_Num_Control","Supporting_Read_Num_Control","Sv_Type"]
        print("\t".join(l_header), file=hout)
        with open(output_prefix + ".savana_sorted.txt", 'r') as hin:
            for line in hin:
                print(line.rstrip('\n'), file=hout)
            
    if not debug:
        os.remove(output_prefix +'.savana_tumor_PASS.bedpe')
        os.remove(output_prefix +'.savana_tumor_repaired.bedpe')
        os.remove(output_prefix +'.savana.vcf.header')
        os.remove(output_prefix +'.savana_filtered.txt')
        os.remove(output_prefix +'.savana_sorted.txt')

    
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

    savanaSVtoBedpe_main(args)
