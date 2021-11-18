#!/usr/bin/python
# Converts a SensitivityProfile file from default format (each base contribution one below the other)
# to a tabular format (one column for each base tuple). 
#
# @author $Author: mario $ 
# @date $Date: 2008-10-01 13:35:53 +0200 (Wed, 01 Oct 2008) $
import sys, csv

reader = csv.reader(open(sys.argv[1], "rb"), delimiter='\t')
# Read all pairs if line does not begin with #
content = [(l, b, c) for l,b,c in reader if not(l.startswith("#"))] 

positionCount = len(set([l for l,b,c in content]))
baseCount = len(content) / positionCount

# Print header
print 'Position',
for base in range(baseCount):
    print content[base * positionCount + 0][1],
print

# Print body
for index in range(positionCount):
    print content[0 + index][0],
    for base in range(baseCount):
        print content[base * positionCount + index][2],
    print
