#!/usr/bin/env python

import sys, subprocess

pdbfn = sys.argv[1]

qt = 0.

for line in open(pdbfn, 'r'):
  if line.startswith('ATOM') or line.startswith('HETATM'):
    sp = line.split()
    q = float(sp[-2])
    qt += q

print 'Total charge:', qt
