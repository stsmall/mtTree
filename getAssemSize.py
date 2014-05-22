#!/usr/bin/env python

import sys

infile = open(sys.argv[1],"r")
lines = infile.readlines()
infile.close()

size = 0

for line in lines:
	if not line.startswith(">"):
		size += len(line.rstrip())

print size
