# vcftobedpeconverter.py
#
# Copyright (c) 2015-2016 Ira Hall lab and The McDonnell Genome Institute
#
# Released under the MIT License.
# https://github.com/hall-lab/svtools/blob/master/LICENSE.txt
# 
# Modified by ken0-1n for python3 compatibility.
#

from .bed import Bed
# from .utils import parse_bnd_alt_string
import re
from builtins import map

class VcfToBedConverter(object):
    '''
    This is a class to take Vcf object(s) and convert them to Bedpe lines
    '''

    def __init__(self):
        '''
        Initialize a new converter
        '''
        pass

    def bnd_breakpoints(self, vcf_variant):
        '''
        Return a tuple containing calculated breakpoints and orientations for a BND variant
        '''
        chrom1 = vcf_variant.chrom
        breakpoint1 = vcf_variant.pos
        orientation1 = '+'

        if vcf_variant.alt.startswith("."):
            orientation1 = '-'
            breakpoint1 -= 1

        return (chrom1,
                breakpoint1,
                breakpoint1,
                orientation1
                )

    @staticmethod
    def adjust_coordinate(vcf_variant, info_tag, start, end):
        '''
        Return adjusted start and end coordinates according to the contents
        of the tag (if it exists)
        '''
        if info_tag in vcf_variant.info:
            span = list(map(int, vcf_variant.info[info_tag].split(',')))
            if len(span) != 2:
                raise ValueError('Invalid value for tag {0}. Require 2 values to adjust coordinates.'.format(info_tag))
            return (start + span[0], end + span[1])
        else:
            return (start, end)

    def convert(self, primary_variant, secondary_variant=None):
        '''
        Convert the passed VCF variant(s) into a BEDPE object
        '''
        vcf_variant = primary_variant

        try:
            sv_type = vcf_variant.info['SVTYPE']
        except KeyError:
            raise ValueError('SVTYPE field required for conversion to BEDPE')

        parser = self.bnd_breakpoints
        c1, s1, e1, o1 = parser(vcf_variant)

        s1, e1 = self.adjust_coordinate(vcf_variant, 'CIPOS', s1, e1)

        orig_name_a = vcf_variant.var_id
        orig_ref_a = vcf_variant.ref
        orig_alt_a = vcf_variant.alt
        info_a = vcf_variant.get_info_string()

        # For MANTA single-ended BNDs, EVENT is not present.
        # XXX This has probably already been calculated outside of this method. May be a candidate to memoize or otherwise cache?
        # By adding to the variant class, perhaps?
        name = vcf_variant.var_id
        if 'EVENT' in vcf_variant.info:
            name = vcf_variant.info['EVENT']
        elif 'MATEID' in vcf_variant.info and vcf_variant.var_id.startswith('Manta'):
            # Specifically handle Manta
            name, end = vcf_variant.var_id.rsplit(':', 1)

        fields = list(map(str, [
            c1,
            max(s1, 0),
            max(e1, 0),
            name,
            vcf_variant.qual,
            o1,
            sv_type,
            vcf_variant.filter,
            orig_name_a,
            orig_ref_a,
            orig_alt_a,
            info_a
            ]))
        if vcf_variant.get_format_string() is not None:
             fields += [vcf_variant.get_format_string(), vcf_variant.get_gt_string()]
        return Bed(fields)
