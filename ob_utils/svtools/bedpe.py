# bedpe.py
#
# Copyright (c) 2015-2016 Ira Hall lab and The McDonnell Genome Institute
#
# Released under the MIT License.
# https://github.com/hall-lab/svtools/blob/master/LICENSE.txt
# 
# Modified by ken0-1n for python3 compatibility.
#

import re
import sys

def find_tag(info_string, tag):
    tag_start = info_string.find(tag)
    while tag_start != -1 and (tag_start > 0 and info_string[tag_start-1] != ';'):
        tag_start = info_string.find(tag, tag_start+1)
    return tag_start


class Bedpe(object):
    def __init__(self, bed_list):
        self.c1 = bed_list[0]
        self.s1 = int(bed_list[1])
        self.e1 = int(bed_list[2])
        self.c2 = bed_list[3]
        self.s2 = int(bed_list[4])
        self.e2 = int(bed_list[5])
        self.name = bed_list[6]
        self.score = self.parse_score(bed_list[7])
        self.o1 = bed_list[8]
        self.o2 = bed_list[9]
        self.svtype = bed_list[10]
        self.filter = bed_list[11]
        self.orig_name1 = bed_list[12]
        self.orig_ref1 = bed_list[13]
        self.orig_alt1 = bed_list[14]
        self.orig_name2 = bed_list[15]
        self.orig_ref2 = bed_list[16]
        self.orig_alt2 = bed_list[17]
        self.malformedFlag = 0
        self.info1 = bed_list[18]
        self.info2 = bed_list[19]
        self.misc = bed_list[20:]
        self.check_malformed()

        # FIXME This is only really needed for varlookup. Something more general would be helpful
        self.cohort_vars = dict()

        try:
            self.svtype = self.retrieve_svtype()
        except ValueError:
            sys.stderr.write('No SVTYPE parseable for {0}'.format('\t'.join(bed_list)))
            sys.exit(1)
        self.af = self.retrieve_af()
        if self.svtype != bed_list[10]:
            sys.stderr.write("SVTYPE at Column 11({0})) and SVTYPE in INFO Column({1}) don't match at variant ID {3}\n".format(str(bed_list[10]), str(self.svtype), self.name))

    @staticmethod
    def parse_score(score):
        if score.isdigit():
            return float(score)
        else:
            return score

    @staticmethod
    def parse_info_tag(info_string, tag):
        '''
        Accessory method to parse out the value of a tag in an info string.
        Make sure to include the equals sign if you are looking for a
        non-boolean tag
        '''

        tag_start = find_tag(info_string, tag)
        if tag_start == -1:
            # If you were looking for a flag then this is the right value.
            # Otherwise your tag doesn't exist. Client code must know how to
            # interpret.
            return False

        tag_end = info_string.find(';', tag_start)
        value_start = tag_start + len(tag)
        if (value_start >= len(info_string)) or (tag_end != -1 and value_start >= tag_end):
            return True
        if tag_end == -1:
            tag_end = None # We didn't find our end index
        return info_string[value_start:tag_end]

    @staticmethod
    def update_info_tag(info_string, tag, new_value):
        '''
        Accessory method to update a tag's value. Like parse_info_tag, make sure to include the equals sign.
        '''

        tag_start = find_tag(info_string, tag)
        if tag_start == -1:
            new_tag = ';' + str(tag);
            new_tag += str(new_value)
            new_info_string = info_string + new_tag
            return new_info_string
            #raise ValueError("Tag {0} doesn't exist".format(tag))

        tag_end = info_string.find(';', tag_start)
        value_start = tag_start + len(tag)
        if (value_start >= tag_end and tag_end != -1) or value_start >= len(info_string):
            raise ValueError("Tag {0} doesn't have a value".format(tag))
        if tag_end == -1:
            tag_end = None # We didn't find our end index
        new_info_string = info_string[:value_start] + new_value
        if tag_end:
            new_info_string += info_string[tag_end:]
        return new_info_string

    @property
    def info(self):
        '''
        Return the appropriate info field if only one is required. Info from the primary variant is preferred if available.
        '''
        if self.info1 == 'MISSING':
            return self.info2
        else:
            return self.info1

    def set_info(self, field, value):
        '''
        Add the info field to the BEDPE line info fields. As BEDPE lines don't know about their headers this is not a safe operation.
        Doesn't add to info field if it is the null character. Probably this is wrong.
        '''
        new_tag = ';' + str(field);
        if value is not None:
            new_tag += '=' + str(value)
        if self.malformedFlag != 1:
            self.info1 = self.info1 + new_tag
        if self.malformedFlag != 2 and self.info2 != '.':
            self.info2 = self.info2 + new_tag

    def check_malformed(self):
        if self.info1 == 'MISSING':
            self.malformedFlag = 1
        if self.info2 == 'MISSING':
            self.malformedFlag = 2

    def retrieve_svtype(self):
        try:
            svtype = re.split('=', ''.join(filter(lambda x: x.startswith('SVTYPE='), self.info.split(';'))))[1]
        except IndexError:
            raise ValueError('SVTYPE field not present in INFO field')
        return svtype

    def retrieve_af(self):
        try:
            af = re.split('=', ''.join(filter(lambda x: x.startswith('AF='), self.info.split(';'))))[1]
        except IndexError:
            af = None
        return af

    def __str__(self):
        '''
        A string representation of the line represented by this object
        '''
        return '\t'.join([
            self.c1,
            str(self.s1),
            str(self.e1),
            self.c2,
            str(self.s2),
            str(self.e2),
            self.name,
            str(self.score),
            self.o1,
            self.o2,
            self.svtype,
            self.filter,
            self.orig_name1,
            self.orig_ref1,
            self.orig_alt1,
            self.orig_name2,
            self.orig_ref2,
            self.orig_alt2,
            self.info1,
            self.info2] +
            self.misc
            )
