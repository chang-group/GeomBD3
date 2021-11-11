#!/usr/bin/python

# Christopher C. Roberts, 2015

import sys


if len(sys.argv) < 2:
  print 'Usage: pdbCenter.py PDBFILENAME'
  sys.exit()

n = 0.
c = [0., 0., 0.]

for line in open(sys.argv[1]):
  line = line.strip()
  if line.startswith('ATOM') or line.startswith('HETATM'):
    c[0] += float(line[30:38])
    c[1] += float(line[38:46])
    c[2] += float(line[46:53])
    n += 1.

for i in range(3):
  c[i] /= n

dx = -c[0]
dy = -c[1]
dz = -c[2]

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

