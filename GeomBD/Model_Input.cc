
#include "Model.h"
#include "Body.h"
#include "Strings.h"
#include "Grid.h"
#include "Grid_ES.h"
#include "Grid_EX.h"
#include "Grid_D.h"
#include "Timer.h"
#include "Session.h"
#include "BindingCriteria.h"
#include "AtomicMass.h"



void Model::parseInputFile() {
  string lfn, line, token, type;
  ifstream cfd;

  cfd.open(ifn.c_str(), ifstream::in);
  while(getline(cfd, line)) {
    if(line[0] == '#') continue;
    parseNextValue(&line, &token);

    if(token.length() > 0) {
      if(token == "grid") {
        parseNextValue(&line, &type);
        parseNextValue(&line, &token);
            if(type == "es") {
              lout << "* Loading electrostatic potential grid \"" << token << "\"...";
              lout.flush();
              Timer *t = new Timer();
              t->start();
              esmaps.push_back(new Grid_ES(token, type));
              t->stop();
              lout << " done. ";
              t->print(&lout);
            } else {
              if(type == "ex") {
                lout << "* Loading exclusion grid \"" << token << "\"...";
                lout.flush();
                Timer *t = new Timer();
                t->start();
                exmaps.push_back(new Grid_EX(token));
                t->stop();
                lout << " done. ";
                t->print(&lout);
              } else {
                if(type == "d") {
                  lout << "* Loading desolvation potential grid \"" << token << "\"...";
                  lout.flush();
                  Timer *t = new Timer();
                  t->start();
                  dmaps.push_back(new Grid_D(token, type));
                  t->stop();
                  lout << " done. ";
                  t->print(&lout);
                } else {
                  lout << "* Loading LJ potential energy grid \"" << token << "\" for atom type \"" << type << "\"...";
                  lout.flush();
                  Timer *t = new Timer();
                  t->start();
                  typemaps.push_back(new Grid_Type(token, type));
                  t->stop();
                  lout << " done. ";
                  t->print(&lout);
                }
              }
            }
      }
      if(token == "debug") {
        parseNextValue(&line, &token);
        lout << "* Loading debug map \"" << token << "\"" << endl;
        debug_map = new Grid_EX(token);
      }
      if(token == "convergence") {
        parseNextValue(&line, &token);
        lout << "* Convergence criteria: " << token << endl;
        convergence = stringToDouble(token);
      }
      if(token == "convWindow") {
        parseNextValue(&line, &token);
        lout << "* Convergence window: " << token << endl;
        convergence_window = stringToInt(token);
      }
      if(token == "threads") {
        parseNextValue(&line, &token);
        __cilkrts_set_param("nworkers", token.c_str());
        lout << "* Attempting to set number of threads to " << token << endl;
      }
      if(token == "temperature") {
        parseNextValue(&line, &token);
        T = stringToDouble(token);
        lout << "* System temperature: " << T << endl;
      }
      if(token == "writeTraj") {
        parseNextValue(&line, &token);
        rate_trj = stringToInt(token);
        lout << "* Writing trajectory every " << rate_trj << " steps." << endl;
      }
      if(token == "writeRestart") {
        parseNextValue(&line, &token);
        rate_rst = stringToInt(token);
        lout << "* Writing restart timers and coordinates every " << rate_rst << " steps." << endl;
      }
      if(token == "useRestart") {
        parseNextValue(&line, &token); // This line should be commented out, right? Copy/paste error?
        if(token=="yes") {
           useRestart = true; 
           lout << "* Simulation starting from restart coordinate (.crd) and timer (.t) files " << bfn << endl;
           }
        else useRestart = false; 
      }
      if(token == "writeBinders") {
        parseNextValue(&line, &bfn);
        lout << "* Writing bound conformations to file " << bfn << endl;
        writeBinders = true;
      }
      if(token == "logBinders") {
        logBinders = true;
      }
      if(token == "logExiters") {
        logExiters = true;
      }
//      if(token == "logDirect") {
//        logDirect = true;
//        boundary = 3; //New 4/16
//      }
      if(token == "betaCalc" or token == "writeLog") {
        parseNextValue(&line, &token);
        rate_beta = stringToInt(token);
        lout << "* Writing association rate information to logfile, every " << rate_beta << " steps." << endl;
      }
      if(token == "timestep") {
        parseNextValue(&line, &token);
        dt_fine = stringToDouble(token);
        parseNextValue(&line, &token);
        dt_coarse = stringToDouble(token);
        parseNextValue(&line, &token);
        dt_scale_start = pow(stringToDouble(token), 2.); //SQUARED VALUE!
        parseNextValue(&line, &token);
        dt_scale_end = pow(stringToDouble(token), 2.);
        lout << "* Timesteps: fine="<<dt_fine<<" ps, coarse=" <<dt_coarse << " ps, scaling from a radius of " << sqrt(dt_scale_start) << "A to a radius of " << sqrt(dt_scale_end) << "A." <<endl;
      }
      if(token == "order") {
        parseNextValue(&line, &token);
        fd_order = stringToInt(token);
        lout << "* Finite difference force approximation using order " << fd_order << "." << endl;
      }
      if(token == "receptor") {
        parseNextValue(&line, &token);
        lout<< "* Receptor filename: " << token << endl;
        if(file_exists(token)) {
          parseReceptorPQR(token);
        } else {
          cout <<"! Receptor file does not exist! Exiting." << endl;
          exit(-1);
        }
      }
      if(token == "fixedReceptor"){
        parseNextValue(&line, &token);
        if(token == "on"){ fixed_receptor = true; }
      }
/*      if(token == "session") {
        parseNextValue(&line, &token);
       
        if(token == "indirect") {
          SessionNAM *sr = new SessionNAM(this);
          sessions.push_back(sr);
          sr->id = sessions.size();
          lout << "* Defining session: indirect contribution to interenzyme transfer" << endl;
        }
      }
*/      if(token == "ligand") {
        SessionA *sr = new SessionA(this); //Automatically create the session; now only one session type in GBD3
        sessions.push_back(sr);
        sr->id = sessions.size();

        parseNextValue(&line, &token); //ligand filename
        lout<< "  + Ligand filename: " << token << endl;
        if(file_exists(token)) {
          parseLigandPQR(token);
        } else {
          cout <<"! Ligand file does not exist! Exiting." << endl;
          exit(-1);
        }
        parseNextValue(&line, &token); //Nreplicates
        sessions[sessions.size()-1]->Nreplicates = stringToInt(token);
        lout<< "  + Ligand replicates: " << token << endl;
      }
      if(token == "ligandStart") {
        parseNextValue(&line, &token);
        if(token == "random") { rand_start = 1; lout << " + Assigning random ligand starting positions" << endl; }
        else {
          string token2, token3, token4;
          parseNextValue(&line, &token2);
          parseNextValue(&line, &token3);
          parseNextValue(&line, &token4);

           if (!token2.empty() && !token3.empty() && !token4.empty()){ //Starts direct/indirect transfer sim if 4 numbers given
             logDirect=true;
             point_start=true;
             pnt_start.x = stringToDouble(token);
             pnt_start.y = stringToDouble(token2);
             pnt_start.z = stringToDouble(token3);
             dircut = stringToDouble(token4);
             cout << "dircut= " << dircut << endl;
             }             
           else if (!token2.empty() && !token3.empty()) { //Starts from an xyz point if 3 numbers given 
             point_start=true;
             pnt_start.x = stringToDouble(token);
             pnt_start.y = stringToDouble(token2);
             pnt_start.z = stringToDouble(token3);
             lout << " + Ligands will start from point " << pnt_start.x <<" "<< pnt_start.y <<" "<< pnt_start.z << endl;
             }
           else if (token2.empty() && token3.empty()) { //Ligands start from a plane at z = bPlnHgt if 1 number given
             usePlane = true;
             bPlnHgt = stringToDouble(token);
             lout << " + Ligands will start from plane at z = " << bPlnHgt << endl;   
             }
           else {  //NAM session is started if 2 numbers given (b sphere start, q sphere exit)
             boundary = 3;
             isNAM = true;
             SessionA *sr = dynamic_cast< SessionA* >(sessions[sessions.size()-1]);
             if(sr) {
               sr->b = stringToDouble(token);
               lout << "  + Starting radius: " << sr->b << " A" << endl;
               sr->q = stringToDouble(token2);
               sr->q2 = sr->q * sr->q;
               lout << "  + Exit radius: " << sr->q << " A" << endl;
               }
             }
           }
        }
      
/*      if(token == "b") {
        parseNextValue(&line, &token);
        SessionNAM *sr = dynamic_cast< SessionNAM* >(sessions[sessions.size()-1]);
        if(sr) {
          sr->b = stringToDouble(token);
          lout << "  + Starting radius: " << sr->b << " A" << endl;
        }
      }
      if(token == "q") {
        parseNextValue(&line, &token);
        SessionNAM *sr = dynamic_cast< SessionNAM* >(sessions[sessions.size()-1]);
        if(sr) {
          sr->q = stringToDouble(token);
          sr->q2 = sr->q * sr->q;
          lout << "  + Exit radius: " << sr->q << " A" << endl;
        }
      }
*/    if(token == "ligandFlux") {
        parseNextValue(&line, &token);
        flux = stringToDouble(token);
        if(flux != 0) lout << " + Applying " << flux << " kcal/mol/A biasing force to ligands in x dimension" << endl;
        }
      if(token == "bindGroups"){
        BindingCriteria *bc = new BindingCriteria(false, true);
        lout << "  + Binding criteria (AND + OR): " << endl;
        int grpID = 0;
        lout << "Group ID is currently: " << grpID << endl; 
        while(parseNextValue(&line, &token)) {
          lout << "Made it into the while loop\n";
          if(token == "OR"){
              grpID++;
              parseNextValue(&line, &token);
              lout << ")\nOR\n(\n";
          } else {
              lout << "AND\n";
          }
          
          lout << "made it past if/else statement\n";
          double bx = stringToDouble(token);
          parseNextValue(&line, &token);
          double by = stringToDouble(token);
          parseNextValue(&line, &token);
          double bz = stringToDouble(token);
          parseNextValue(&line, &token);
          int laid = stringToInt(token);
          parseNextValue(&line, &token);
          lout << "I stand before the mighty r assignment\n";
          double r = stringToDouble(token);
          bc->addPair(bx, by, bz, laid-1, r, grpID); //addPair function adds bind criteria into vector 'bindingSet' within vector 'pairs'
          lout << "    - LAID: " << laid << " within " << r << "A of (" << bx << " " << by << " " << bz << ")" << endl;
          if(laid > sessions[sessions.size()-1]->conformations[0]->beads.size()) { //New
             cout << " ! Error: LAID out of range" << endl; 
             exit(0); 
             }
        }
        sessions[sessions.size()-1]->bindingCriteria.push_back(bc); //bindingCriteria size is always 1 since outside of while loop here
        lout << bc->pairs.size() << " binding criteria (bindingSets) found\n";
        if(!useRestart) {
           for(int i=0; i<bc->pairs.size(); i++) { bc->critNbind.push_back(0.); bc->critTimes.push_back(0.); }
           }  
      }
      if(token == "ceiling"){
      parseNextValue(&line, &token);
      ceiling = stringToDouble(token);                      
      }
      if(token == "randomStart"){
        parseNextValue(&line, &token);
        if(token == "on"){ rand_start = 1; }
      }  
      if(token == "boundary"){
        parseNextValue(&line, &token);
        if(token == "rect") { boundary = 1; lout << "Using rectangular prism boundaries\n"; }
        else if(token == "hex") { boundary = 2; lout << "Using hexagonal prism boundaries\n"; }
        else if(token == "sphere") { 
           if (usePlane) { cout << " ! Error: planar starting conditions must be used with rectangular or hexagonal boundaries\n"; std::exit(0); }
           else { boundary = 4; lout << "Using spherical boudary conditions\n"; }
           }
        if (isNAM or logDirect) { boundary = 3; lout << "Using spherical NAM boundary conditions\n"; } //Last check forces NAM q sphere if isNAM or transfer simulation
        }            
      if(token == "boundWall"){
        parseNextValue(&line, &token);
        wall = stringToDouble(token);
       // if(logDirect) { //New way to force q exit for direct/indirect by assign
         // SessionA *sr = dynamic_cast< SessionA* >(sessions[sessions.size()-1]);
           // if(sr) {
             // sr->q = wall;
              //sr->q2 = sr->q * sr->q;
              //lout << "  + Exit radius: " << sr->q << " A" << endl;
              //}
          //}
      }
      if(token == "planeWidth"){
        parseNextValue(&line, &token);
        planeWidth = stringToDouble(token);
      }
//      if(token == "directCut"){  //New 4/15
 //       parseNextValue(&line, &token);
//        dircut = stringToDouble(token);
//      }
      if(token == "maxSims") {
        parseNextValue(&line, &token);
        max_simulations = stringToInt(token);
        lout << "* Setting maximum number of completed replicate simulations to " << max_simulations << "." << endl;
      }
    }
  }

  cfd.close();

  if(logDirect) { //New way to force q exit for direct/indirect by assign 4/16
    SessionA *sr = dynamic_cast< SessionA* >(sessions[sessions.size()-1]);
    if(sr) { 
       sr->q = wall; sr->q2 = sr->q * sr->q; lout << "  + Exit radius: " << sr->q << " A" << endl; }
       sr->b = dircut; lout << "  + Indirect radius: " << sr->b << " A" << endl; //New 0513
    }

  //Calculate simulation volume
  if(boundary==1) volume = 4*wall*wall*ceiling;
  else if(boundary==2) { double a = tan(M_PI/6)*2*wall; volume = sqrt(3)*1.5*a*a*ceiling; }
  else if(boundary==4) volume = (4*M_PI*pow(wall,3)) / 3;
  if(volume!=0) lout << "* Simulation volume is " << volume << " A^3" << endl;

  // Determine system geometry based on largest grid 3/2/21
  for(int i=0; i < esmaps.size(); i++) {
    if(esmaps[i]->origin[0] < bounds_min.x) bounds_min.x = esmaps[i]->origin[0];
    if(esmaps[i]->origin[1] < bounds_min.y) bounds_min.y = esmaps[i]->origin[1];
    if(esmaps[i]->origin[2] < bounds_min.z) bounds_min.z = esmaps[i]->origin[2];
    bounds_max.x = esmaps[i]->origin[0] + ( (esmaps[0]->N[0]-1) * esmaps[0]->delta );
    bounds_max.y = esmaps[i]->origin[1] + ( (esmaps[0]->N[1]-1) * esmaps[0]->delta );
    bounds_max.z = esmaps[i]->origin[2] + ( (esmaps[0]->N[2]-1) * esmaps[0]->delta );
  }
  for(int i=0; i < typemaps.size(); i++) {
    if(typemaps[i]->origin[0] < bounds_min.x) bounds_min.x = typemaps[i]->origin[0];
    if(typemaps[i]->origin[1] < bounds_min.y) bounds_min.y = typemaps[i]->origin[1];
    if(typemaps[i]->origin[2] < bounds_min.z) bounds_min.z = typemaps[i]->origin[2];
    if(typemaps[i]->origin[0] + ( (typemaps[0]->N[0]-1) * typemaps[0]->delta ) > bounds_max.x) {
       bounds_max.x = typemaps[i]->origin[0] + ( (typemaps[0]->N[0]-1) * typemaps[0]->delta );
       }
    if(typemaps[i]->origin[1] + ( (typemaps[0]->N[1]-1) * typemaps[0]->delta ) > bounds_max.y) {
       bounds_max.y = typemaps[i]->origin[1] + ( (typemaps[0]->N[1]-1) * typemaps[0]->delta );
       }
    if(typemaps[i]->origin[2] + ( (typemaps[0]->N[2]-1) * typemaps[0]->delta ) > bounds_max.z) {
       bounds_max.z = typemaps[i]->origin[2] + ( (typemaps[0]->N[2]-1) * typemaps[0]->delta );
       }
  }
  for(int i=0; i < exmaps.size(); i++) {
    if(exmaps[i]->origin[0] < bounds_min.x) bounds_min.x = exmaps[i]->origin[0];
    if(exmaps[i]->origin[1] < bounds_min.y) bounds_min.y = exmaps[i]->origin[1];
    if(exmaps[i]->origin[2] < bounds_min.z) bounds_min.z = exmaps[i]->origin[2];
    if(exmaps[i]->origin[0] + ( (exmaps[0]->N[0]-1) * exmaps[0]->delta ) > bounds_max.x) {
       bounds_max.x = exmaps[i]->origin[0] + ( (exmaps[0]->N[0]-1) * exmaps[0]->delta );
       }
    if(exmaps[i]->origin[1] + ( (exmaps[0]->N[1]-1) * exmaps[0]->delta ) > bounds_max.y) {
       bounds_max.y = exmaps[i]->origin[1] + ( (exmaps[0]->N[1]-1) * exmaps[0]->delta );
       }
    if(exmaps[i]->origin[2] + ( (exmaps[0]->N[2]-1) * exmaps[0]->delta ) > bounds_max.z) {
       bounds_max.z = exmaps[i]->origin[2] + ( (exmaps[0]->N[2]-1) * exmaps[0]->delta );
       }
  }

  receptor_max.x = exmaps[0]->origin[0] + ( (exmaps[0]->N[0]-1) * exmaps[0]->delta ) ;
  receptor_min.x = exmaps[0]->origin[0] ;
  receptor_max.y = exmaps[0]->origin[1] + ( (exmaps[0]->N[1]-1) * exmaps[0]->delta ) ;
  receptor_min.y = exmaps[0]->origin[1] ;
  receptor_max.z = exmaps[0]->origin[2] + ( (exmaps[0]->N[2]-1) * exmaps[0]->delta ) ;
  receptor_min.z = exmaps[0]->origin[2] ;
  
// BEGIN new output of bounds min/max
  lout << "Grid min/max extent = " << bounds_min.x << " " << bounds_min.y <<" " << bounds_min.z << " "
       << bounds_max.x << " " << bounds_max.y<< " " << bounds_max.z << endl;
  lout << "* Receptor min/max extent:"  << receptor_min.x << " " << receptor_min.y << " " << receptor_min.z << " / "
       << receptor_max.x << " " << receptor_max.y << " " << receptor_max.z << endl;
  lout << "Receptor center at: " << center.x <<" "<< center.y <<" "<< center.z << endl;
// END new output
} //Close parseInputFile()

void Model::parseReceptorPQR(string rfn) {
  string line, token;
  vector<double> rx;
  vector<double> ry;
  vector<double> rz;
  vector<double> radii;
  double cr[3] = { 0., 0., 0. };

  ifstream fd(rfn);

  while(getline(fd, line)) {
    if(starts_with(&line, "ATOM")) {
      string at_raw = line.substr(12, 4);
      string at = trim(at_raw);
      double rd = stringToDouble(line.substr(69, 6));
      radii.push_back(rd);
      double x = stringToDouble(line.substr(30, 10));
      double y = stringToDouble(line.substr(40, 10));
      double z = stringToDouble(line.substr(50, 10));      
      rx.push_back(x);
      ry.push_back(y);
      rz.push_back(z);

      x_pen_chk.push_back(x); // New
      y_pen_chk.push_back(y);
      z_pen_chk.push_back(z);

      cr[0] += x;
      cr[1] += y;
      cr[2] += z;
    }
  }
  cr[0] /= rx.size();
  cr[1] /= rx.size();
  cr[2] /= rx.size();

  if(rx.size() == 1) {
    receptor_radius = radii[0];
  } else {
    double rmsd = 0., dr = 0.;

    for(int i=0; i < rx.size(); i++) {
      dr = rx[i] - cr[0];
      rmsd += dr * dr;
      dr = ry[i] - cr[1];
      rmsd += dr * dr;
      dr = rz[i] - cr[2];
      rmsd += dr * dr;
    }

    rmsd /= rx.size();
    receptor_radius = sqrt(rmsd);
  }


  center.x = cr[0];
  center.y = cr[1];
  center.z = cr[2];

  rx.clear();
  ry.clear();
  rz.clear();
}


void Model::parseLigandPQR(string lfn) {
  string line, token;
  Body *bi = new Body(this, sessions[sessions.size()-1]);
  Bead *bj = NULL;

  int Nconfs = 0;

  ifstream fd(lfn);

  while(getline(fd, line)) {
    if(starts_with(&line, "ATOM") or starts_with(&line, "HETATM")) {
      if(bi == NULL) {
        bi = new Body(this, sessions[sessions.size()-1]);
      }
      bj = new Bead();

      char element = line[13];
      double x = stringToDouble(line.substr(30, 10));
      double y = stringToDouble(line.substr(40, 10));
      double z = stringToDouble(line.substr(50, 10));
      double q = stringToDouble(line.substr(60, 8));    // Here is where ligand q is read. Changing (70, 7) to other numbers TC 5/31
      string at_raw = line.substr(12, 4);
      string at = trim(at_raw);
      bj->r = stringToDouble(line.substr(69, 6));
      bj->m = atomicMass(at);

      if(bj->m == 0.) {
        lout << "! FATAL: No mass for ligand bead type: " << at << endl;
        exit(-1);
      }

      bj->R.x = x;
      bj->R.y = y;
      bj->R.z = z;
      bj->q = q;
      bj->type = at;

      //New solpar assignment
      if (at=="H") bj->S = 0.00051 + 0.01097*fabs(bj->q);
      else if (at=="C") bj->S = -0.00143 + 0.01097*fabs(bj->q);
      else if (at=="N") bj->S = -0.00162 + 0.01097*fabs(bj->q);
      else if (at=="O") bj->S = -0.00251 + 0.01097*fabs(bj->q);
      else if (at=="S") bj->S = -0.00214 + 0.01097*fabs(bj->q);
      else if (at=="P") bj->S = -0.00110 + 0.01097*fabs(bj->q);
      else bj->S = -0.00110 + 0.01097*fabs(bj->q);
      bj->S *= 0.1; //Desolvation re-weighting coef
      bi->beads.push_back(bj);
    }
    if(starts_with(&line, "END")) {
      if(bi != NULL) {
        sessions[sessions.size()-1]->conformations.push_back(bi);
        bi->define();
        bi = NULL;
        Nconfs++;
      }
    }
  }
  lout << "    - Loaded " << Nconfs << " ligand conformation" << endl;
}






