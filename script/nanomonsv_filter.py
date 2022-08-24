#! /usr/bin/env python
import sys, csv

def basic_filter(row):
    filter_flag = False
    if row["Chr_1"].endswith("decoy"): filter_flag = True
    if row["Chr_2"].endswith("decoy"): filter_flag = True
    if row["Supporting_Read_Num_Control"] != "0": filter_flag = True
    if row["Chr_1"] == row["Chr_2"] and row["Dir_1"] == '+' and row["Dir_2"] == '-':
        sv_size = int(row["Pos_2"]) - int(row["Pos_1"]) + len(row["Inserted_Seq"]) - 1
        if sv_size < 100: filter_flag = True
        # if row["Is_Simple_Repeat"] != "None" and row["Insert_Type"] in ["None", "---"]: filter_flag = True
        # insert_type = row["Insert_Type"]
        # if insert_type not in ["None", "---"] and row["Is_Simple_Repeat"] != "None": filter_flag = True
    return filter_flag

input_file = sys.argv[1]

svkey2read_num = {}
with open(input_file, 'r') as hin:
    for row in csv.DictReader(hin, delimiter = '\t'):
        if basic_filter(row): continue
        svkey = '\t'.join([row[x] for x in ["Chr_1", "Pos_1", "Dir_1", "Chr_2", "Pos_2", "Dir_2"]])
        if row["Inserted_Seq"] == "---": row["Inserted_Seq"] = ""
        svkey2read_num[svkey] = (int(row["Supporting_Read_Num_Control"]), len(row["Inserted_Seq"]))

header_print = False
svkey2writtern = {}
with open(input_file, 'r') as hin:
    for row in csv.DictReader(hin, delimiter = '\t'):
        if header_print == False:
            print('\t'.join(row))
            header_print = True
        if basic_filter(row): continue
        inseq_len_self = 0 if row["Inserted_Seq"] == "---" else len(row["Inserted_Seq"])
        svkey_line_self = '\t'.join([row[x] for x in ["Chr_1", "Pos_1", "Dir_1", "Chr_2", "Pos_2", "Dir_2"]])
        if svkey_line_self in svkey2writtern: continue
        match_flag1 = False
        for svkey_line in svkey2read_num:
            if svkey_line == svkey_line_self: continue
            svkey = svkey_line.split('\t')
            if row["Chr_1"] == svkey[0] and row["Chr_2"] == svkey[3]:
                if row["Dir_1"] == svkey[2] and row["Dir_2"] == svkey[5]:
                    if abs(int(row["Pos_1"]) - int(svkey[1])) <= 20 and abs(int(row["Pos_2"]) - int(svkey[4])) <= 20:
                        if int(row["Supporting_Read_Num_Control"]) < svkey2read_num[svkey_line][0]: 
                            match_flag1 = True
                        elif int(row["Supporting_Read_Num_Control"]) == svkey2read_num[svkey_line][0]:
                            if inseq_len_self > svkey2read_num[svkey_line][1]:
                                match_flag1 = True
                            elif inseq_len_self == svkey2read_num[svkey_line][1]:
                                if svkey_line < svkey_line_self:
                                    match_flag1 = True 
        match_flag2 = False
        if match_flag1 == False and match_flag2 == False:       
            print('\t'.join(row.values()))
            svkey2writtern[svkey_line_self] = 1

