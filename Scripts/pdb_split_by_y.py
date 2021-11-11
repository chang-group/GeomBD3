#!/usr/bin/python

# Christopher C. Roberts, 2015

import sys

ycutoff = float(sys.argv[2])
abovef = open(sys.argv[3],'w')
belowf = open(sys.argv[4],'w')

for line in open(sys.argv[1]):
  if line.startswith('ATOM') or line.startswith('HETATM'):
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    if y > ycutoff:
      abovef.write(line)
    else:
      belowf.write(line)

