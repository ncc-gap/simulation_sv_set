import sys
REFERENCE = sys.argv[1]
START = sys.argv[2]
END = sys.argv[3]

print_flag = False
with open(REFERENCE, 'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>' + START):
            print_flag = True
        elif line.startswith('>' + END):
            break
    
        if print_flag:
            print(line)
