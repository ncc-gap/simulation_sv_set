import sys

in_file = sys.argv[1]
out_prefix  = sys.argv[2]
with open(in_file,'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>'):
            line = line+"_"+out_prefix
        print(line)
