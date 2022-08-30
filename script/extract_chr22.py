import sys
REFERENCE = sys.argv[1]

print_flag = False
with open(REFERENCE, 'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>chr22'):
            print_flag = True
        elif line.startswith('>chrX'):
            break
    
        if print_flag:
            print(line)
