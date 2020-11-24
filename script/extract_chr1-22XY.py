import sys

in_file = sys.argv[1]

print_flag = False
with open(in_file,'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>chr1'):
            print_flag = True
        elif line.startswith('>chrM'):
            break
    
        if print_flag:
            print(line)
