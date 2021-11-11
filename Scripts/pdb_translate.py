#!/usr/bin/python

# Christopher C. Roberts, 2015

import sys

if len(sys.argv) < 5:
  print 'Usage: pdbTranslate.py PDBFILENAME DX DY DZ'
  sys.exit()

dx = float(sys.argv[2])
dy = float(sys.argv[3])
dz = float(sys.argv[4])

print 'REMARK translating', dx, dy, dz

for line in open(sys.argv[1]):
  line = line.strip()
  if line.startswith('ATOM') or line.startswith('HETATM'):
    x = '%8.3f' % (float(line[30:38]) + dx)
    y = '%8.3f' % (float(line[38:46]) + dy)
    z = '%8.3f' % (float(line[46:54]) + dz)
    newline = line[:30] + x + y +z + line[54:]
    print newline
  else:
    print line

