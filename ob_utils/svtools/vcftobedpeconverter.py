# vcftobedpeconverter.py
#
# Copyright (c) 2015-2016 Ira Hall lab and The McDonnell Genome Institute
#
# Released under the MIT License.
# https://github.com/hall-lab/svtools/blob/master/LICENSE.txt
# 
# Modified by ken0-1n for python3 compatibility.
#

from .bedpe import Bedpe
# from .utils import parse_bnd_alt_string
import re
from builtins import map

class VcfToBedpeConverter(object):
    '''
    This is a class to take Vcf object(s) and convert them to Bedpe lines
    '''

    def __init__(self, strands_tag, end_tag):
        '''
        Initialize a new converter
        '''
        #pass
        self.STRANDS_TAG = strands_tag
        self.END_TAG = end_tag

    def bnd_breakpoints(self, vcf_variant):
        '''
        Return a tuple containing calculated breakpoints and orientations for a BND variant
        '''
        chrom1 = vcf_variant.chrom
        breakpoint1 = vcf_variant.pos
        orientation1 = orientation2 = '+'
        sep, chrom2, breakpoint2 = parse_bnd_alt_string(vcf_variant.alt)
        breakpoint2 = int(breakpoint2)

        if vcf_variant.alt.startswith(sep):
            orientation1 = '-'
            breakpoint1 -= 1

        if sep == '[':
            orientation2 = '-'
            breakpoint2 -= 1

        return (chrom1,
                breakpoint1,
                breakpoint1,
                chrom2,
                breakpoint2,
                breakpoint2,
                orientation1,
                orientation2)

    #@staticmethod
    #def simple_breakpoints(vcf_variant):
    def simple_breakpoints(self, vcf_variant):
        '''
        Return a tuple containing breakpoints and orientations for simple SVs
        '''
        breakpoint1 = vcf_variant.pos
        try:
            breakpoint2 = int(float(vcf_variant.info[self.END_TAG]))
        except KeyError:
            raise ValueError('END entry in VCF required for conversion to BEDPE')

        orientation1 = '+'
        orientation2 = '-'

        #if 'STRANDS' in vcf_variant.info:
        #    strands = vcf_variant.info['STRANDS']
        if self.STRANDS_TAG in vcf_variant.info:
            strands = vcf_variant.info[self.STRANDS_TAG]
            orientation1, orientation2 = strands[:2]

        return (vcf_variant.chrom,
            breakpoint1,
            breakpoint1,
            vcf_variant.chrom,
            breakpoint2,
            breakpoint2,
            orientation1,
            orientation2)

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
        if primary_variant is None:
            vcf_variant = secondary_variant

        try:
            sv_type = vcf_variant.info['SVTYPE']
        except KeyError:
            raise ValueError('SVTYPE field required for conversion to BEDPE')

        parser = self.simple_breakpoints
        if sv_type == 'BND':
            parser = self.bnd_breakpoints

        c1, s1, e1, c2, s2, e2, o1, o2 = parser(vcf_variant)

        s1, e1 = self.adjust_coordinate(vcf_variant, 'CIPOS', s1, e1)
        s2, e2 = self.adjust_coordinate(vcf_variant, 'CIEND', s2, e2)

        orig_name_a = vcf_variant.var_id
        orig_ref_a = vcf_variant.ref
        orig_alt_a = vcf_variant.alt
        info_a = vcf_variant.get_info_string()
        if primary_variant is None:
            info_a = "MISSING"
            orig_name_a = orig_ref_a = orig_alt_a = '.'
            c1, s1, e1, o1, c2, s2, e2, o2 = c2, s2, e2, o2, c1, s1, e1, o1

        info_b = '.'
        orig_name_b = orig_ref_b = orig_alt_b = '.'
        if sv_type == 'BND':
            if secondary_variant is None:
                info_b = "MISSING"
            else:
                info_b = secondary_variant.get_info_string()
                orig_name_b = secondary_variant.var_id
                orig_ref_b = secondary_variant.ref
                orig_alt_b = secondary_variant.alt
                sc1, ss1, se1, sc2, ss2, se2, so1, so2 = parser(secondary_variant)
                s2, e2 = self.adjust_coordinate(secondary_variant, 'CIPOS', ss1, se1)

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
            c2,
            max(s2, 0),
            max(e2, 0),
            name,
            vcf_variant.qual,
            o1,
            o2,
            sv_type,
            vcf_variant.filter,
            orig_name_a,
            orig_ref_a,
            orig_alt_a,
            orig_name_b,
            orig_ref_b,
            orig_alt_b,
            info_a,
            info_b,
            ]))
        if vcf_variant.get_format_string() is not None:
             fields += [vcf_variant.get_format_string(), vcf_variant.get_gt_string()]
        return Bedpe(fields)

def parse_bnd_alt_string(alt_string):
    '''
    Parse the BND alt string and return separators and region
    '''
    # NOTE The below is ugly but intended to match things like [2:222[ and capture the brackets
    result = re.findall(r'([][])(.+?)([][])', alt_string)
    assert result, "%s\n" % alt_string
    #sys.stderr.write("%s\n" % alt_string)
    sep1, region, sep2 = result[0]
    assert sep1 == sep2
    chrom2, breakpoint2 = region.rsplit(':', 1)
    breakpoint2 = breakpoint2
    return sep1, chrom2, breakpoint2
