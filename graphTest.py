#path = "OutputTest.txt"
##with open(path,'r') as f:
##    for line in f:
##        print(line)
#lines = [line for line in open(path)]
#print(lines[1])
#print(lines[0])
#print(lines[1].strip().split(","))

import csv
#print(dir(csv))
path = "OutputTest.txt"
file = open(path,newline='')        # To work across platforms (or so I hear)
readfunc = csv.reader(file)
header = next(readfunc)             # First line in header
data = [row for row in readfunc]    # Read remaining data

print(header)
print(data[0])