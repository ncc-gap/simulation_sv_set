#! /usr/bin/env python
import sys, pysam, gzip, math
import os

input_file = sys.argv[1]
golden_data_bedpe = sys.argv[2]
check_margin = 200 # int(sys.argv[3])
insertion_margin = 0.2 # float(sys.argv[4])
output = sys.argv[3]

golden_data_db = pysam.TabixFile(golden_data_bedpe)


l_del = []
hout = open(output, 'w')
with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        tchr1, tpos1, tdir1, tchr2, tpos2, tdir2 = F[0], F[1], F[2], F[3], F[4], F[5]
        insert_seq = F[6]
        insert_size = 0 if insert_seq == "---" else len(insert_seq)

        golden_data_flag = False
        tabix_error_flag = False
        try:
            records = golden_data_db.fetch(F[0], max(0, int(tpos1) - 200), int(tpos1) + 200)
        except Exception as e:
            #print(e)
            tabix_error_flag = True

        if not tabix_error_flag:
            for record_line in records:
                record = record_line.split('\t')
                # nanomonsv = ins, golden data = dup
                if tchr1 == record[0] and tchr2 == record[3] and \
                ((tdir1 == "+" and  tdir2 ==  "-"  and record[8] == "-" and record[9] == "+") or \
                (tdir1 == "*" and tdir2 == "*")):

                    if int(tpos1) >= int(record[1]) - check_margin and \
                    int(tpos1) <= int(record[2]) + check_margin and \
                    int(tpos2) >= int(record[4]) - check_margin - insert_size + 1 and \
                    int(tpos2) <= int(record[5]) + check_margin - insert_size + 1 : 

                        golden_data_flag = True
                        break
                    
                    if int(tpos1) >= int(record[1]) - check_margin + insert_size  and \
                    int(tpos1) <= int(record[2]) + check_margin + insert_size and \
                    int(tpos2) >= int(record[4]) - check_margin and \
                    int(tpos2) <= int(record[5]) + check_margin : 
                        golden_data_flag = True
                        break
                    
                # nanomonsv = ins, golden data = ins
                if tchr1 == record[0] and tchr2 == record[3] and \
                ((tdir1 == record[8] and tdir2 == record[9]) or (tdir1 == "*" and tdir2 == "*")) and \
                record[10] == "INS":
                
                    record_insert_size = int(record[11])
                    ins_diff = record_insert_size - insert_size
                    ins_diff_per = (abs(ins_diff) / record_insert_size)
                    if ins_diff_per <= insertion_margin:

                        if int(tpos1) >= int(record[1]) - check_margin and \
                        int(tpos1) <= int(record[2]) + check_margin and \
                        int(tpos2) >= int(record[4]) - check_margin and \
                        int(tpos2) <= int(record[5]) + check_margin: 
                            golden_data_flag = True
                            break
                    
                if tchr1 == record[0] and tchr2 == record[3] and \
                (tdir1 == record[8] and tdir2 == record[9]):
                    
                    insert_size1 = insert_size
                    insert_size2 = insert_size
                    if tdir1 == "+": insert_size1 = insert_size*-1
                    if tdir1 == "-": insert_size1 = insert_size
                    if tdir2 == "+": insert_size2 = insert_size*-1
                    if tdir2 == "-": insert_size2 = insert_size

                    if int(tpos1) >= int(record[1]) - check_margin and \
                    int(tpos1) <= int(record[2]) + check_margin and \
                    int(tpos2) >= int(record[4]) - check_margin + insert_size2 and \
                    int(tpos2) <= int(record[5]) + check_margin + insert_size2: 
                        golden_data_flag = True
                        break
                        
                    elif int(tpos1) >= int(record[1]) - check_margin + insert_size1 and \
                    int(tpos1) <= int(record[2]) + check_margin + insert_size1 and \
                    int(tpos2) >= int(record[4]) - check_margin  and \
                    int(tpos2) <= int(record[5]) + check_margin: 
                        golden_data_flag = True
                        break

                if tchr1 == record[0] and tchr2 == record[3] and \
                (tdir1 == "*" and tdir2 == "*"):
                    
                    if int(tpos1) >= int(record[1]) - check_margin and \
                    int(tpos1) <= int(record[2]) + check_margin and \
                    int(tpos2) >= int(record[4]) - check_margin - insert_size and \
                    int(tpos2) <= int(record[5]) + check_margin + insert_size: 
                        golden_data_flag = True
                        break
                        
                    elif int(tpos1) >= int(record[1]) - check_margin - insert_size and \
                    int(tpos1) <= int(record[2]) + check_margin + insert_size and \
                    int(tpos2) >= int(record[4]) - check_margin  and \
                    int(tpos2) <= int(record[5]) + check_margin: 
                        golden_data_flag = True
                        break

        print_line = "\t".join(F)
        if golden_data_flag:
            l_record = [record[0],record[2],record[8],record[3],record[5],record[9],record[10],record[11]]
            print_line = print_line +"\t"+"\t".join(l_record)
            
            l_record_key = [record[0],record[1],record[2],record[3],record[4],record[5],record[8],record[9]]
            l_del.append(",".join(l_record_key))
        else:
            l_empty = ["" for _ in range(7)]
            print_line = print_line +"\t"+"\t".join(l_empty)
        print(print_line, file=hout)
hout.close()

output_prefix, ext = os.path.splitext(output)
hout = open(output_prefix + ".false_negative.txt", "w")
with gzip.open(golden_data_bedpe, "rt") as hin:
    for line in hin:
        line = line.rstrip('\n')
        F = line.split('\t')
        key = ",".join(F[:6])+','+F[8]+','+F[9]
        if key not in l_del:
            print('\t'.join(F), file=hout)
hout.close()


