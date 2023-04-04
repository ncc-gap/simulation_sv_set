import sys
REFERENCE = sys.argv[1]

print_flag = False
with open(REFERENCE, 'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>chr1'):
            print_flag = True
        elif line.startswith('>chrM'):
            break
    
        if print_flag:
            print(line)
