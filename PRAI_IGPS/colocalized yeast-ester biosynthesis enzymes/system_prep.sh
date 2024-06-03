#!/bin/bash

export GBD3=/path_of_directory/GBD3_master  #change the path of GBD3_master

#Parameterize ligand (crp) and receptor (IGPS) PQR files from PDB files
python2 $GBD3/bin/Parameterize.py -d $GBD3/Parameters/AA.gdbp -i crp.pdb -o crp.pqr
python2 $GBD3/bin/Parameterize.py -d $GBD3/Parameters/AA.gdbp -i IGPS.pdb -o IGPS.pqr

#Run Gridder programs to create electrostatic, Lennard-Jones, ligand desolvation, and volume exclusion grids for the system
$GBD3/bin/Gridder-ES -d $GBD3/Parameters/AA.gdbp -r IGPS.pqr -o IGPS- -n 10
$GBD3/bin/Gridder-LJ -d $GBD3/Parameters/AA_0.2.gdbp -r IGPS.pqr -l crp.pqr -o IGPS- -n 10 -p 12
$GBD3/bin/Gridder-D -d $GBD3/Parameters/AA.gdbp -r IGPS.pqr -o IGPS- -n 10 -p 12
$GBD3/bin/Gridder-EX -r IGPS.pqr -o IGPS- -n 10

