#!/usr/bin/python
# For every pair "A\tB" in file DICTIONARY, replaces A with B in INPUTFILE
#
# Mario Fasold, last update 2009/06/08
import sys

## If wrong number of command line arguments supplied, print error message and exit
if len(sys.argv) < 2:
    print "Usage: %s DICTIONARY [INPUTFILE]"%sys.argv[0]
    sys.exit()

## If no INPUTFILE given, use STDIN
inputfile = sys.stdin
if len(sys.argv) >= 3:
    inputfile = open(sys.argv[2])

def repleaceAll(s, replaces):
    for match, repl in replaces:
        s = s.replace(match, repl)
    return s

# Read key-value pairs into a dictionary
replaces = {}
for line in open(sys.argv[1], "rb"):
    if line[0] == '#': # omit comments
        continue
    lineSplit = line.strip().split('\t')
    if len(lineSplit) > 1:
        replaces[lineSplit[0]] = lineSplit[1]
    else: # replace empty/undefined variables as well
        replaces[lineSplit[0]] = "" 

# Convert to tuple-list and sort by length(string-to-replace) to avoid earlier replace by substrings!
replaceTuples = replaces.items()
replaceTuples.sort(lambda x,y:cmp(len(x[0]), len(y[0]) ), reverse=1)

# Replace all the terms in each line
for line in inputfile:
   print repleaceAll(line, replaceTuples).rstrip()



