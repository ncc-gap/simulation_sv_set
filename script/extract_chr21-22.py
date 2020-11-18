

print_flag = False
with open('GRCh38.d1.vd1.fa','r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>chr21'):
            print_flag = True
        elif line.startswith('>chrX'):
            break
    
        if print_flag:
            print(line)
