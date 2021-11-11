#!/usr/bin/python

# Christopher C. Roberts, 2015

import sys
from math import cos, sin


if len(sys.argv) < 5:
  print 'Usage: pdbRotate.py PDBFILENAME DTHETAX DTHETAY DTHETAZ'
  sys.exit()



dx = float(sys.argv[2])
dy = float(sys.argv[3])
dz = float(sys.argv[4])



def rotate(v, dx, dy, dz):
  rm = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
  cost = 0.
  sint = 0.

  p = [0., 0., 0.]

  if dx != 0.:
    cost = cos(dx)
    sint = sin(dx)
    rm[0][0] = 1
    rm[0][1] = 0.
    rm[0][2] = 0.

    rm[1][0] = 0.
    rm[1][1] = cost
    rm[1][2] = -sint

    rm[2][0] = 0.
    rm[2][1] = sint
    rm[2][2] = cost

    for n in range(3):
      for m in range(3):
        p[n] += rm[n][m] * v[m]

    for n in range(3):
      v[n] = p[n]
      p[n] = 0.

  if dy != 0.:
    cost = cos(dy)
    sint = sin(dy)
    rm[0][0] = cost
    rm[0][1] = 0.
    rm[0][2] = sint

    rm[1][0] = 0.
    rm[1][1] = 1.
    rm[1][2] = 0.

    rm[2][0] = -sint
    rm[2][1] = 0.
    rm[2][2] = cost

    for n in range(3):
      for m in range(3):
        p[n] += rm[n][m] * v[m]

    for n in range(3):
      v[n] = p[n]
      p[n] = 0.

  if dz != 0.:
    cost = cos(dz)
    sint = sin(dz)
    rm[0][0] = cost
    rm[0][1] = -sint
    rm[0][2] = 0.

    rm[1][0] = sint
    rm[1][1] = cost
    rm[1][2] = 0.

    rm[2][0] = 0.
    rm[2][1] = 0.
    rm[2][2] = 1.

    for n in range(3):
      for m in range(3):
        p[n] += rm[n][m] * v[m]

    for n in range(3):
      v[n] = p[n]
      p[n] = 0.

  return v


for line in open(sys.argv[1]):
  line = line.strip()
  if line.startswith('ATOM') or line.startswith('HETATM'):
    v = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    v = rotate(v, dx, dy, dz)
    x = '%8.3f' % (v[0])
    y = '%8.3f' % (v[1])
    z = '%8.3f' % (v[2])
    newline = line[:30] + x + y +z + line[54:]
    print newline

