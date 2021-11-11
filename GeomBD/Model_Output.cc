#include "Model.h"
#include "Strings.h"
#include "BindingCriteria.h"
#include "Session.h"

void writeBead(Bead *bi, unt index, char chain, fstream &outf) {
  outf << "ATOM  ";
  outf.width(5);
  outf << right << index;
  outf << "      ";
  outf.width(3);
  outf << " P " << " ";
  outf << chain;
  outf << "        ";
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.x/1.;
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.y/1.;
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.z/1.;
  outf << " ";
  outf.setf(ios::right, ios::adjustfield);
  outf.precision(4);
  outf.width(7);
  outf << bi->q << " " << bi->type << endl;
}



void Model::writeCoordinatesPQR() {
  string ofnpqr = ofn;
  if (step==0) ofnpqr.append(".pqr"); //cout << "Successfully wrote starting coordinates .pqr file\n";}
  else ofnpqr.append(".crd"); // New 4/13 erased brackets
  // create output file
  fstream outf(ofnpqr.c_str(), ios::out/*| ios::app*/);

  if(!outf) //EDIT
    cout << "The files was NOT opened succesfully\n"; //EDIT
  if(outf) //EDIT
    //cout << "goodbit is set for output file\n"; //EDIT

  int index = 0;
  char chain = 'A';

  for(int b=0; b < ligands.size(); b++) {
    Body *ligand = ligands[b];
  //  cout << "Created pointer to Body object in ligands[for_counter] vector\n"; //EDIT
    ligand->writePDB(outf, chain);
    chain++;
    if(chain > 'z') chain = 'A';
  }
  outf << "END" << endl;
  outf.close();
}

// Testing new timer outpout 12/2
void Model::writeTimes() {
  string tofn = ofn;
  tofn.append(".t");
  fstream outf(tofn.c_str(), ios::out);
   for (int t=0; t<ligands.size(); t++) {
      Body *ligand = ligands[t];
      outf << ligand->t << endl;
      }

  BindingCriteria *bcout = sessions[0]->bindingCriteria[0];
   for(int i=0; i<bcout->critNbind.size(); i++) {
      outf << bcout->critNbind[i] << endl; 
      outf << bcout->critTimes[i] << endl;
      }
   for(int i=0; i<sessions[0]->bindtimes.size(); i++) {
      outf << sessions[0]->bindtimes[i] << endl;
      }
  outf.close();
}

void Model::openTrajectoryDCD() {
  //cout << "About to call writeCoordinatesPQR function\n"; //EDIT
  writeCoordinatesPQR(); //need a "topology" for DCD loading
  //cout << "writeCoordinatesPQR called succesfully\n"; //EDIT

  for(int b=0; b < ligands.size(); b++) {
    Body *ligand = ligands[b];
    Natoms += ligand->N;
  }

  int it = 0;
  float ft = 0.;
  char buffer[256];
  string dofn = ofn; //New 3/10

  dofn.append(".dcd"); //New 3/10
  ofd.open(dofn.c_str()/*ofn.c_str()*/, ios::out | ios::binary);
  if(! ofd.is_open()) {
    cout << "! Can't open DCD trajectory for writing." << endl;
    return;
  }

  //cout << "About to write some sort of header to the DCD file\n"; //EDIT

  ofd.seekp(0, std::ios::end);
  /*  0*/it = 84;                   ofd.write((char*)&it, sizeof(it));
  /*  4*/strcpy(buffer, "CORD");    ofd.write(buffer, 4);
  /*  8*/it = 1;/*nframe*/          ofd.write((char*)&it, sizeof(it));
  /* 12*/it = 1;                    ofd.write((char*)&it, sizeof(it));
  /* 16*/it = 1;                    ofd.write((char*)&it, sizeof(it));
  /* 20*/it = 1;/*nframe*/          ofd.write((char*)&it, sizeof(it));
  /* 24*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 28*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 32*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 36*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 40*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 44*/ft = 0.;                   ofd.write((char*)&ft, sizeof(ft));
  /* 48*/it = 1;                    ofd.write((char*)&it, sizeof(it));
  /* 52*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 56*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 60*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 64*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 68*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 72*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 76*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 80*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 84*/it = 24;                   ofd.write((char*)&it, sizeof(it));
  /* 88*/it = 84;                   ofd.write((char*)&it, sizeof(it));
  /* 92*/it = 164;                  ofd.write((char*)&it, sizeof(it));
  /* 96*/it = 2;                    ofd.write((char*)&it, sizeof(it));
  /*100*/sprintf(buffer, "%160s\0", "REMARKS"); ofd.write(buffer, 160);
  /*260*/it = 164;                  ofd.write((char*)&it, sizeof(it));
  /*264*/it = 4;                    ofd.write((char*)&it, sizeof(it));
  /*268*/it = Natoms;               ofd.write((char*)&it, sizeof(it));
  /*272*/it = 4;                    ofd.write((char*)&it, sizeof(it));
}


void Model::writeCoordinatesDCD() {
  double dt;
  float ft;
  int it;
  char buffer[256];

  /*  0*/it = 48;                                    ofd.write((char*)&it, sizeof(it));
  /*  4*/dt = 1000.;                                 ofd.write((char*)&dt, sizeof(dt));
  /* 12*/strcpy(buffer, "        ");                 ofd.write(buffer, 8);
  /* 20*/dt = 1000.;                                 ofd.write((char*)&dt, sizeof(dt));
  /* 28*/strcpy(buffer, "                ");         ofd.write(buffer, 16);
  /* 44*/dt = 1000.;                                 ofd.write((char*)&dt, sizeof(dt));
  /* 52*/it = 48;                                    ofd.write((char*)&it, sizeof(it));
  for(int c=0; c < 3; c++) {
    /* 56*/it = Natoms * 4;                          ofd.write((char*)&it, sizeof(it));
    for(int l=0; l < ligands.size(); l++) {
      Body *ligand = ligands[l];
      for(int b=0; b < ligand->beads.size(); b++) {
        if(c == 0) ft = ligand->beads[b]->R.x;
        if(c == 1) ft = ligand->beads[b]->R.y;
        if(c == 2) ft = ligand->beads[b]->R.z;
        /* 60*/ ofd.write((char*)&ft, sizeof(ft));
      }
    }
    /*  0*/it = Natoms * 4;                          ofd.write((char*)&it, sizeof(it));
  }

  Nframes++;
}


void Model::closeTrajectoryDCD() {
  float ft;
  int it;
  char buffer[256];

  ofd.seekp(0, std::ios::end);
  /*  0*/it = 84;                   ofd.write((char*)&it, sizeof(it));
  /*  4*/strcpy(buffer, "CORD");    ofd.write(buffer, 4);
  /*  8*/it = Nframes;              ofd.write((char*)&it, sizeof(it));
  /* 12*/it = 1;                    ofd.write((char*)&it, sizeof(it));
  /* 16*/it = 1;                    ofd.write((char*)&it, sizeof(it));
  /* 20*/it = Nframes;              ofd.write((char*)&it, sizeof(it));
  /* 24*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 28*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 32*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 36*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 40*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 44*/ft = 0.;                   ofd.write((char*)&ft, sizeof(ft));
  /* 48*/it = 1;                    ofd.write((char*)&it, sizeof(it));
  /* 52*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 56*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 60*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 64*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 68*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 72*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 76*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 80*/it = 0;                    ofd.write((char*)&it, sizeof(it));
  /* 84*/it = 24;                   ofd.write((char*)&it, sizeof(it));
  /* 88*/it = 84;                   ofd.write((char*)&it, sizeof(it));
  /* 92*/it = 164;                  ofd.write((char*)&it, sizeof(it));
  /* 96*/it = 2;                    ofd.write((char*)&it, sizeof(it));
  /*100*/sprintf(buffer, "%160s\0", "REMARKS"); ofd.write(buffer, 160);
  /*260*/it = 164;                  ofd.write((char*)&it, sizeof(it));
  /*264*/it = 4;                    ofd.write((char*)&it, sizeof(it));
  /*268*/it = Natoms;               ofd.write((char*)&it, sizeof(it));
  /*272*/it = 4;                    ofd.write((char*)&it, sizeof(it));

  ofd.close();
}
