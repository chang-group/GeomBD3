#!/usr/bin/python

import sys, math

maxc = [0., 0., 0.]
minc = [1e90, 1e90, 1e90]
coords = []
center = [0., 0., 0.]
N = 0


def calc():
  center[0] /= N
  center[1] /= N
  center[2] /= N

  #invSumRij = 0.
  #for i in range(N):
  #  for j in range(i+1, N):
  #    drx = coords[i][0] - coords[j][0]
  #    dry = coords[i][1] - coords[j][1]
  #    drz = coords[i][2] - coords[j][2]
  #    invSumRij += 1. / math.sqrt(drx*drx + dry*dry + drz*drz)
  #invRhyd = (1.0 / pow(N, 2)) * invSumRij
  #Rhyd = 1.0 / invRhyd

  rmsd = 0.

  for c in coords:
    for i in range(3):
      dr = c[i] - center[i]
      rmsd += dr * dr
      maxc[i] = max([c[i], maxc[i]])
      minc[i] = min([c[i], minc[i]])

  rmsd /= N
  rmsd = math.sqrt(rmsd)
  print 'ROJ:', rmsd
  print 'Center:', center
  print 'Max:', maxc
  print 'Min:', minc
  #print 'Rhyd:', Rhyd

pqr = True
if sys.argv[1][-3:].lower() == 'pdb':
  pqr = False

for line in open(sys.argv[1], 'r'):
  if line.startswith('ATOM'):
    if pqr:
      x = float(line[30:40])
      y = float(line[40:50])
      z = float(line[50:60])
    else:
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
    coords.append([x,y,z])
    center[0] += x
    center[1] += y
    center[2] += z
    N += 1
  if line.startswith('END'):
    calc()
    maxc = [0., 0., 0.]
    minc = [1e90, 1e90, 1e90]
    coords = []
    center = [0., 0., 0.]
    N = 0


