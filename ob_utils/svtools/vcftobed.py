# vcftobedpe.py
#
# Copyright (c) 2015-2016 Ira Hall lab and The McDonnell Genome Institute
#
# Released under the MIT License.
# https://github.com/hall-lab/svtools/blob/master/LICENSE.txt
#
# Modified by ken0-1n for python3 compatibility.
# 

import argparse
import sys
import time

from .vcf.file import Vcf 
from .vcf.variant import Variant
from .vcftobedconverter import VcfToBedConverter

def vcfToBed(vcf_file, bed_out):
    converter = VcfToBedConverter()
    vcf = Vcf()
    in_header = True
    header = []
    sample_list = []
    v = []
    for line in vcf_file:
        if in_header:
            if line[0:2] == '##':
                if line.split('=')[0] == '##fileformat':
                    line = '##fileformat=' + "BEDPE" + '\n'
                if line.split('=')[0] == '##fileDate':
                    line = '##fileDate=' + time.strftime('%Y%m%d') + '\n'
                header.append(line)
                continue
            elif line[0] == '#' and line[1] != '#':
                sample_list = line.rstrip().split('\t')[9:]
                header.append(line)
                continue
            else:
                # print header
                in_header = False
                vcf.add_header(header)
                if "SVTYPE" in [info.id for info in vcf.info_list]:
                   vcf.add_info_after("SVTYPE", "POS", 1, 'Integer', 'Position of the variant described in this record')
                header=vcf.get_header()
                bed_out.write(header[:header.rfind('\n')] + '\n')
                final_header_line = ['#CHROM_A',
                        'START_A',
                        'END_A',
                        'ID',
                        'QUAL',
                        'STRAND_A',
                        'TYPE',
                        'FILTER',
                        'NAME_A',
                        'REF_A',
                        'ALT_A',
                        'INFO_A']

                if len(sample_list) > 0:
                    bed_out.write('\t'.join(final_header_line + ['FORMAT','\t'.join(map(str,sample_list))]) + '\n')
                else:
                    bed_out.write('\t'.join(final_header_line) + '\n')

        v = line.rstrip().split('\t')
        var = Variant(v, vcf)
        var.set_info("POS", var.pos)
        # If there is no MATEID then assume this is a single-ended BND and simply output
        if var.info['SVTYPE'] != 'BND' or 'MATEID' not in var.info:
            bed_out.write(str(converter.convert(var)) + '\n')
        else:
            sys.stderr.write('Warning: ' + var.var_id + ' is not the single end SV \n')

    # close the files
    bed_out.close()
    return

def run_vcf2bed(input_vcf, output_bedpe):
    
    with open(output_bedpe, 'w') as hout:
        if input_vcf.endswith('.gz'):
            import gzip
            with gzip.open(input_vcf, 'rb') as hin:
                return vcfToBed(hin, hout)
        else:
            with open(input_vcf, 'r') as hin:
                return vcfToBed(hin, hout)

