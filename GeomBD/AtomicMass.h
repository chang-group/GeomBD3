#ifndef _atomic_mass
#define _atomic_mass

inline double atomicMass(string element) {
  if(element == "H") return 1.0079;
  if(element == "He") return 4.0026;
  if(element == "Li") return 6.941;
  if(element == "Be") return 9.0122;
  if(element == "B") return 10.811;
  if(element == "C") return 12.0107;
  if(element == "N") return 14.0067;
  if(element == "O") return 15.9994;
  if(element == "F") return 18.9984;
  if(element == "Ne") return 20.1797;
  if(element == "Na") return 22.9897;
  if(element == "Mg") return 24.305;
  if(element == "Al") return 26.9815;
  if(element == "Si") return 28.0855;
  if(element == "P") return 30.9738;
  if(element == "S") return 32.065;
  if(element == "Cl") return 35.453;
  if(element == "K") return 39.0983;
  if(element == "Ar") return 39.948;
  if(element == "Ca") return 40.078;
  if(element == "Sc") return 44.9559;
  if(element == "Ti") return 47.867;
  if(element == "V") return 50.9415;
  if(element == "Cr") return 51.9961;
  if(element == "Mn") return 54.938;
  if(element == "Fe") return 55.845;
  if(element == "Ni") return 58.6934;
  if(element == "Co") return 58.9332;
  if(element == "Cu") return 63.546;
  if(element == "Zn") return 65.39;
  if(element == "Ga") return 69.723;
  if(element == "Ge") return 72.64;
  if(element == "As") return 74.9216;
  if(element == "Se") return 78.96;
  if(element == "Br") return 79.904;
  if(element == "Kr") return 83.8;
  if(element == "Rb") return 85.4678;
  if(element == "Sr") return 87.62;
  if(element == "Y") return 88.9059;
  if(element == "Zr") return 91.224;
  if(element == "Nb") return 92.9064;
  if(element == "Mo") return 95.94;
  if(element == "Tc") return 98;
  if(element == "Ru") return 101.07;
  if(element == "Rh") return 102.9055;
  if(element == "Pd") return 106.42;
  if(element == "Ag") return 107.8682;
  if(element == "Cd") return 112.411;
  if(element == "In") return 114.818;
  if(element == "Sn") return 118.71;
  if(element == "Sb") return 121.76;
  if(element == "I") return 126.9045;
  if(element == "Te") return 127.6;
  if(element == "Xe") return 131.293;
  if(element == "Cs") return 132.9055;
  if(element == "Ba") return 137.327;
  if(element == "La") return 138.9055;
  if(element == "Ce") return 140.116;
  if(element == "Pr") return 140.9077;
  if(element == "Nd") return 144.24;
  if(element == "Pm") return 145;
  if(element == "Sm") return 150.36;
  if(element == "Eu") return 151.964;
  if(element == "Gd") return 157.25;
  if(element == "Tb") return 158.9253;
  if(element == "Dy") return 162.5;
  if(element == "Ho") return 164.9303;
  if(element == "Er") return 167.259;
  if(element == "Tm") return 168.9342;
  if(element == "Yb") return 173.04;
  if(element == "Lu") return 174.967;
  if(element == "Hf") return 178.49;
  if(element == "Ta") return 180.9479;
  if(element == "W") return 183.84;
  if(element == "Re") return 186.207;
  if(element == "Os") return 190.23;
  if(element == "Ir") return 192.217;
  if(element == "Pt") return 195.078;
  if(element == "Au") return 196.9665;
  if(element == "Hg") return 200.59;
  if(element == "Tl") return 204.3833;
  if(element == "Pb") return 207.2;
  if(element == "Bi") return 208.9804;
  if(element == "Nu") return 12.00; //New test ligand
  if(element == "Po") return 12.00; //New test ligand
  return 0;
}

#endif
