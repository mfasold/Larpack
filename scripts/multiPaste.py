#!/usr/bin/python
# Works alike the linux command 'paste' but, for each line, 
# the command "multipaste j i file1 file2 file3 .. fileX
#  * prints the j-th column of file1 
#  * print's i-th columns of files 1...X
#
# This version can handle more files, but has a large memory footprint
#
# What if nothing works: first "paste" all files. then "cut" every second column.
import sys
import csv
selectedLeadColumn = int(sys.argv[1]) - 1
selectedColumn = int(sys.argv[2]) - 1

filenames = sys.argv[3:]

# Create array of vectors in first run
new_table = list()
for row in csv.reader(open(filenames[1],"rb"), delimiter='\t'):
    new_table.append(["x" for i in range(len(filenames) + 1)])
    new_table[-1][0] = row[0]

# Populate array
for file_index, filename in enumerate(filenames):
    for row_index, row in enumerate(csv.reader(open(filename, "rb"), delimiter='\t')):
        new_table[row_index][file_index+1] = row[1] 

for line in new_table:
    print "\t".join(line)


