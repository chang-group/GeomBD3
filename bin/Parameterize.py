#!/usr/bin/env python2

import sys, subprocess



def at_synonym_key(at):
  if at == 'H5\'': return 'H5A'
  elif at == 'H5\'\'': return 'H5B'
  elif at == 'H2\'': return 'H2A'
  elif at == 'H2\'\'': return 'H2B'
  elif at == 'HO3\'': return 'H3T'
  elif at == 'HO5\'': return 'H5T' 
  elif at == 'HO2\'': return 'HO2A'
  elif at == 'HO2\'\'': return 'HO2B'
  return at
def at_synonym(at):
  if at == 'H5\'': return 'H5\'1'
  elif at == 'H5\'\'': return 'H5\'2'
  elif at == 'H2\'': return 'H2\'1'
  elif at == 'H2\'\'': return 'H2\'2'
  elif at == 'HO5\'': return 'HO5\'1'
  elif at == 'HO5\'\'': return 'HO5\'2'
  elif at == 'HO2\'': return 'HO2\'1'
  elif at == 'HO2\'\'': return 'HO2\'2'
  return at


def at_to_element(at):
  return at[0]

def process_receptor(outf, af, ff, resname, resid, resdata):
  debug = open('debug', 'a')  #EDIT
  debug.write("we are inside proccess receptor. The current resname is ")
  debug.write(str(resname))
  debug.write("\n")
  if resname in ff.keys():
    debug.write("WE are inside the first if statement\n")
    for line in resdata:
      at = line[12:16].strip()
      at = at_synonym_key(at)
      el = at_to_element(at)
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
      debug.write(str(el))
      try:
        param = ff[resname][at]
        #print '%s  1.00  0.00   %8.4f %s' % (line[:54], param[0], param[1])
        #TODO - I'm sure this is where the error in writing out to file and filling hex values with blanks is
        pref = '%s%s%s' % (line[:12], str(' %s' % param[1].ljust(3, ' ')), line[16:30])
        outf.write('%s%10.4f%10.4f%10.4f %7.4f %6.4f\n' % (pref, x, y, z, param[0], af[param[1]]['r'])) #PQR
      except KeyError:
        print 'REMARK No charge assign for atom type', at
        pref = '%s%s%s' % (line[:12], str(' %s' % el.ljust(3, ' ')), line[16:30])
        #print '%s  1.00  0.00   %8.4f %s' % (line[:54], 0, at[0])
        outf.write('%s%10.4f%10.4f%10.4f %7.4f %6.4f\n' % (pref, x, y, z, 0, af[at[0]]['r'])) #PQR
  else:
    debug.write("We didnt make the if sttement, so we are in else\n")
    outf.write('REMARK No forcefield assignment for the following residue. Using OpenBABEL\n')
    tmpfd = open('/tmp/het.pdb', 'w')
    for line in resdata:
      tmpfd.write(line)
    tmpfd.close()
    subprocess.call(["babel", "-ipdb", "/tmp/het.pdb", "-opdb", "/tmp/heth.pdb"])
    subprocess.call(["babel", "-ipdb", "/tmp/heth.pdb", "-omol2", "/tmp/het.mol2"])
    Q = []
    tmpfd = open('/tmp/het.mol2', 'r')
    tmpfd.readline()
    tmpfd.readline()
    n = int(tmpfd.readline().split()[0])
    line = tmpfd.readline()
    while not line.startswith("@<TRIPOS>ATOM"):
      line = tmpfd.readline()
    for line in tmpfd:
      if line[0] == '@' or line.strip() == '': break
      sp = line.split()
      q = float(sp[-1])
      Q.append(q)
    i = 0
    debug.write("we made it to the last for loop\n")
    for line in open('/tmp/heth.pdb', 'r'):
      if line.startswith('ATOM') or line.startswith('HETATM'):
        at = at_synonym(line[12:16].strip())
        debug.write("the element symbol read in is: ")
        debug.write(str(at))
        debug.write("\n")
        el = at_to_element(at)
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        debug.write("we made it past variable assignments\n")
        pref = '%s%s%s' % (line[:12], str(' %s' % el.ljust(3, ' ')), line[16:30])
        debug.write("We made it past the pref string creation\n")
        debug.write("pref     x     y     z     Q[i]     el     af     af[el][r]\n")
        debug.write(str(pref))
        debug.write(",")
        debug.write(str(x))
        debug.write(",")
        debug.write(str(y))
        debug.write(",")
        debug.write(str(z))
        debug.write(",")
        debug.write(str(Q[i]))
        debug.write(",")
        debug.write(str(el))
        debug.write(",") 
        #debug.write(str(af))
        #debug.write(",")
        debug.write(str(af[el]))
        debug.write(",")
        debug.write(str(af[el]['r']))
        debug.write("\n")
        outf.write('%s%10.4f%10.4f%10.4f %7.4f %6.4f\n' % (pref, x, y, z, Q[i], af[el]['r'])) #PQR
        debug.write("we made it past the write to file\n")
        i += 1
  debug.write("We made it to the end of the else statement\n")
  debug.close() #EDIT



def run(paramfn, pdbfn, outfn):
  ff = {}
  af = {}
  debug = open('debug', 'w')  #EDIT
  debug.write(" ")
  debug.close()

  for line in open(paramfn, 'r'):
    sp = line.split()
    if len(sp) < 5: continue
    if sp[0] == 'frag_par':
      try:
        ff[sp[1]][sp[2]] = [float(sp[3]), sp[4]]
      except KeyError:
        ff[sp[1]] = {}
        ff[sp[1]][sp[2]] = [float(sp[3]), sp[4]]
        debug = open('debug','a')
        debug.write("We never make it past run into proccess receptor for gold\n")
        debug.write(str(sp[1]))
        debug.write("\n")
        debug.close()
    if sp[0] == 'atom_par':
      debug = open('debug','a')
      debug.write(str(sp[1]))
      debug.write("\n")
      debug.close()
      af[sp[1]] = {'r':float(sp[2])/2., 'eps':float(sp[3]), 'vol':float(sp[4]), 'sol':float(sp[5])} #LOOK HERE
      debug = open('debug','a')
      debug.write("we make it past the exception thrown and go for another loop\n")
      debug.close()

  resn = None
  resi = -6969
  resd = []
  term = False

  outf = open(outfn, 'w')

  for line in open(pdbfn, 'r'):
    if line.startswith('ATOM') or line.startswith('HETATM'):
      #EDIT - added in nested elif statements to catch hex value March 15, 2018 2121hrs
      if 'a' in line[23:26]:
        rid = int(line[23:26],16)
      elif 'b' in line[23:26]:
        rid = int(line[23:26],16)
      elif 'c' in line[23:26]:
        rid = int(line[23:26],16)
      elif 'd' in line[23:26]:
        rid = int(line[23:26],16)
      elif 'e' in line[23:26]:
        rid = int(line[23:26],16)
      elif 'f' in line[23:26]:
        rid = int(line[23:26],16)
      else:
        rid = int(line[23:26]) 
      #END-EDIT
      #rid = int(line[23:26])
      rnm = line[17:20].strip()
      if rid != resi:
        if resi == -6969:
          tnm = 'N%s' % rnm
          if tnm in ff.keys(): rnm = tnm
        else:
          process_receptor(outf, af, ff, resn, resi, resd)
          resd = []
          if term:
            outf.write('TER\n')
        if term:
          tnm = 'N%s' % rnm
          if tnm in ff.keys(): rnm = tnm
          term = False
        resi = rid
        resn = rnm
      resd.append(line)
      if line[12:16].strip() == 'OXT':
        resn = 'C%s' % resn
        term = True
    if line.startswith('END'):
      if len(resd) != 0:
        process_receptor(outf, af, ff, resn, resi, resd)
        resi = -6969
        resd = []
        resn = None
        outf.write('END\n')


if __name__ == '__main__':
  if len(sys.argv) < 7:
    print 'Usage:', sys.argv[0], '-d [Parm.gbdp]', '-i [Molecule.PDB]', '-o [Molecule.PQR]'
  if len(sys.argv) == 7:
    argd = { sys.argv[1]: sys.argv[2], sys.argv[3]: sys.argv[4], sys.argv[5]: sys.argv[6] }
    run(argd['-d'], argd['-i'], argd['-o'])
