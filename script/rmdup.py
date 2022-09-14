#! /usr/bin/env python

import sys, csv

input_file = sys.argv[1]
from collections import namedtuple

def get_sv_class(F):

    tclass = None
    tinseq = '' if F["Inserted_Seq"] == '---' else F["Inserted_Seq"]
    if F["Chr_1"] != F["Chr_2"]:
        tclass = "Translocation"
    elif F["Dir_1"] == F["Dir_2"]:
        tclass = "Inversion"
    elif F["Dir_1"] == '-' and F["Dir_2"] == '+':
        tclass = "Duplication" 
    elif F["Dir_1"] == '+' and F["Dir_2"] == '-':
        if len(tinseq) > int(F["Pos_2"]) - int(F["Pos_1"]) + 1:
            tclass = "Insertion"
        else:
            tclass = "Deletion"

    return(tclass)


def basic_filter(row):

    filter_flag = False
    if row["Chr_1"].endswith("decoy"): filter_flag = True
    if row["Chr_2"].endswith("decoy"): filter_flag = True
    #if row["Supporting_Read_Num_Control"] != "0": filter_flag = True
    if not row["Supporting_Read_Num_Control"] in ["---", "."] and row["Supporting_Read_Num_Control"] != "0": filter_flag = True

    if row["Chr_1"] == row["Chr_2"] and row["Dir_1"] == '+' and row["Dir_2"] == '-':
        sv_size = int(row["Pos_2"]) - int(row["Pos_1"]) + len(row["Inserted_Seq"]) - 1
        if sv_size < 100: filter_flag = True


    return filter_flag

def get_checked_read_num_tumor(F):
        checked_read_num_tumor = "" 
        if F["Checked_Read_Num_Tumor"] == "---" or F["Checked_Read_Num_Tumor"] == "." :
            checked_read_num_tumor = "1" 
        else:
            checked_read_num_tumor = F["Checked_Read_Num_Tumor"]
        return checked_read_num_tumor

svkey2info = {}
with open(input_file, 'r') as hin:
    for F in csv.DictReader(hin, delimiter = '\t'):

        tclass = get_sv_class(F)
        tinseq = '' if F["Inserted_Seq"] == '---' else F["Inserted_Seq"]

        svkey = '\t'.join([F[x] for x in ["Chr_1", "Pos_1", "Dir_1", \
            "Chr_2", "Pos_2", "Dir_2"]])
        svkey2info[svkey] = [int(get_checked_read_num_tumor(F)), len(tinseq), tclass]

svkey2writtern = {}
with open(input_file, 'r') as hin:
    dreader = csv.DictReader(hin, delimiter = '\t')
    print('\t'.join(dreader.fieldnames))
    for row in dreader:

        if basic_filter(row): continue
        tclass = get_sv_class(row) 

        inseq_len_self = 0 if row["Inserted_Seq"] == "---" else len(row["Inserted_Seq"])
        svkey_line_self = '\t'.join([row[x] for x in ["Chr_1", "Pos_1", "Dir_1", "Chr_2", "Pos_2", "Dir_2"]])
        if svkey_line_self in svkey2writtern: continue

        match_flag1 = False
        for svkey_line in svkey2info:
            svkey = svkey_line.split('\t')
            if row["Chr_1"] == svkey[0] and row["Chr_2"] == svkey[3]:
                if row["Dir_1"] == svkey[2] and row["Dir_2"] == svkey[5]:
                    if abs(int(row["Pos_1"]) - int(svkey[1])) <= 50 and abs(int(row["Pos_2"]) - int(svkey[4])) <= 50:

                        if int(get_checked_read_num_tumor(row)) < svkey2info[svkey_line][0]:
                            match_flag1 = True
                        elif int(get_checked_read_num_tumor(row)) == svkey2info[svkey_line][0]:
                            if inseq_len_self > svkey2info[svkey_line][1]:
                                match_flag1 = True
                            elif inseq_len_self == svkey2info[svkey_line][1]:
                                if svkey_line > svkey_line_self:
                                    match_flag1 = True

            # check co-existence of insertion and duplication
            if tclass == "Insertion" and row["Chr_1"] == svkey[0] and row["Chr_2"] == svkey[3]:
                if svkey[2] == '-' and svkey[5] == '+': 
                    if abs(int(row["Pos_1"]) - int(svkey[1])) <= 50 or abs(int(row["Pos_2"]) - int(svkey[4])) <= 50:
                        if 0.8 * (int(svkey[4]) - int(svkey[1])) <= inseq_len_self <=  1.2 * (int(svkey[4]) - int(svkey[1])):
                            match_flag1 = True

            # if row["Pos_1"] == "15708165" and svkey[1] == "15707904":
            #     import pdb; pdb.set_trace()

            # check duplicated insertion sites
            if tclass == "Insertion" and row["Chr_1"] == svkey[0] and row["Chr_2"] == svkey[3]:
                if svkey2info[svkey_line][2] == "Insertion" and abs(int(row["Pos_1"]) - int(svkey[1])) <= 1000:
                    if 0.8 * svkey2info[svkey_line][1] <= inseq_len_self <= 1.2 * svkey2info[svkey_line][1]:
                        if int(row["Pos_1"]) > int(svkey[1]): match_flag1 = True

        if match_flag1 == False:
            print('\t'.join(row.values()))
            svkey2writtern[svkey_line_self] = 1



