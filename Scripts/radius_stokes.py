#!/usr/bin/env python
import math, sys

#D = 0.062
D = float(sys.argv[1])
T = 298.
pi = math.pi
v = 2.39e-25
Na = 6.0221415e23
kb = 1.9858775e-3 / Na

print 'Rhyd =', ((kb * T) / (6 * pi * v * D)), 'A'
