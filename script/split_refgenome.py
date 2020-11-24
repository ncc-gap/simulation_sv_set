
import sys 
in_fa = sys.argv[1]
chromosome = sys.argv[2]
suffix = sys.argv[3]

print_flag = False
with open(in_fa,'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>'+chromosome+suffix):
            print_flag = True
        elif print_flag and line.startswith('>'):
            break
    
        if print_flag:
            print(line)
