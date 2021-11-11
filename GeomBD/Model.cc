#include "Session.h"
#include "Model.h"
#include "Timer.h"
#include "BindingCriteria.h"
#include "Body.h"


Model::Model() {
  writeBinders = false;

  logBinders = false;
  logExiters = false;

  rngCPU = NULL;
  rand = NULL;

  T = 298.;
  viscosity = 1.0/*cP*/ * 0.001/*cP->Pa.s->J.s/m^3*/ * 2.39e-4 /*J->kcal*/ * 1e12 /*s->ps*/ * 1e-30 /*1/m^3->1/A^3*/ * Na /*kcal->kcal/mol*/;
  threads = 1;
  fd_order = 2;
  rate_trj = 10000;
  rate_beta = 5000;
  convergence = 1e-4;
  convergence_window = 100;
  max_simulations = -1;

  debug_map = NULL;

  center.x = center.y = center.z = 0.;
  system_extent = 0.;
  receptor_radius = 1e20;
  bounds_min.x = bounds_min.y = bounds_min.z = 1e9;
  bounds_max.x = bounds_max.y = bounds_max.z = -1e9;

  step = 0;
  done = false;
  dt_fine = 0.010;
  dt_coarse = 1.000;
  dt_scale_start = 100.;
  dt_scale_end = 500.;

  Natoms = 0;
  Nframes = 0;
}


Model::Model(string inputfn, string outputfn, string logfn, string rstcrd, string rsttime) : Model() {
//Model::Model(string inputfn, string outputfn, string logfn) : Model() {
  ifn = inputfn;
  if(!file_exists(ifn)) {
    cout << "! Input file \"" << ifn << "\" does not exist. Exiting." << endl;
    exit(EXIT_FAILURE);
  }
  
  ofn = outputfn;
  lfn = logfn;
  lfn.append(".log"); // New 3/10
  lout.open(lfn, ios_base::out);

  parseInputFile();  
  threads = __cilkrts_get_nworkers();

  //cout << "Thread count is: " << threads << endl; //EDIT
  //cout << "We managed to parse the input file\n"; //EDIT
  
  rsc = rstcrd; //New
  rst = rsttime; //New

  //if(file_exists(rsc) or file_exists(rst)) useRestart = true;

  if(useRestart) {
     if(!file_exists(rsc) or !file_exists(rst)) {
        cout << "! Error: Simulation is set as a restart, but one or more restart files are not specified\n";
        exit(EXIT_FAILURE);
        }
     }
  if(!useRestart) {
     if(file_exists(rsc) or file_exists(rst)) {
        cout << "! Error: Simulation is set as a non-restart simulation, but one or more restart files have been specified\n"
             << "! If this is a restart, set useRestart = yes\n";
        exit(EXIT_FAILURE);
        }
     }

// Reading in restart files
  if (useRestart) {
    ifstream irst;
    string irstline;
    irst.open(rst, ios::in);
    int ci=0;
    while (irst) {
       getline (irst, irstline);
       if(!irst) {break;}
       else { 
         if(ci < sessions[sessions.size()-1]->Nreplicates)  vrt.push_back(atof(irstline.c_str()));
         if(ci >= sessions[sessions.size()-1]->Nreplicates) vrb.push_back(atof(irstline.c_str())); 
         ci++;
         }
      }
    irst.close();
    cout << "vrt size = " << vrt.size() << " vrb size = " << vrb.size() << " Contents of vectors vrt and vrb:\n";
    for (int i=0; i<vrt.size(); i++) { cout << vrt[i] << endl; }
    for (int i=0; i<vrb.size(); i++) { cout << vrb[i] << endl; }
    // Assigning restart binding info to bindingCriteria
    BindingCriteria *bco = sessions[sessions.size()-1]->bindingCriteria[0];
    for (int i=0; i<bco->pairs.size(); i++) { 
       bco->critNbind.push_back(vrb[i*2]); 
       cout << "pushed back " << vrb[i*2] << endl;
       cout << "value in critNbind: " << bco->critNbind[i] << endl;
       bco->critTimes.push_back(vrb[i*2+1]);
       cout << "pushed back " << vrb[i*2+1] << endl;
       cout << "value in critTimes: " << bco->critTimes[i] << endl;
       }
    for (int i=bco->pairs.size()*2; i<vrb.size(); i++) {
       sessions[0]->bindtimes.push_back(vrb[i]);
       cout << "restart bindtimes " << i-bco->pairs.size()*2 << " = " << sessions[0]->bindtimes[i-bco->pairs.size()*2] << endl;
       }
    cout << "  bindtimes vec size= " << sessions[0]->bindtimes.size() << endl;
    ifstream irsc;
    string irscline;
    irsc.open(rsc, ios::in);
    while (irsc) {
       getline (irsc, irscline);
       if(!irsc) {break;}
       else {
         if (starts_with(&irscline, "ATOM")) { 
           vrc.push_back(atof(irscline.substr(30, 8).c_str()));
           vrc.push_back(atof(irscline.substr(38, 8).c_str()));
           vrc.push_back(atof(irscline.substr(46, 8).c_str()));
           }
         }
       }
  irsc.close();
  } // New end of if useRestart read rst files

  populateLigands();

  initializeRNG(time(NULL));

  //cout << "initializeRNG was called succesfully\n"; //EDIT
}



Model::~Model() {
  cout << "The destructor is called, delete has yet to be called\n"; //EDIT
  for(int k=0; k < threads; k++) 
    vslDeleteStream(&rngCPU[k]);
  cout << "Deleteion succesful inside the destructo\n"; //EDIT
}



void Model::initializeRNG(int seed) {
  lout << "* Initializing random number generator (N=" << threads << ")" << endl;

  srand(seed);

  rngCPU = (VSLStreamStatePtr*)calloc(threads, sizeof(VSLStreamStatePtr));
  for(int k=0; k < threads; k++) {
    vslNewStream(&rngCPU[k], VSL_BRNG_MT19937, seed+k);
  }

  rand = (vertex*)calloc(2 * ligands.size(), sizeof(vertex));
}



void Model::generateNormal() {
  int blockSize = ((2*ligands.size()*3) / threads);
  cilk_for(int k=0; k < threads; k++) {
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, rngCPU[k], blockSize, (double*)rand+(k*blockSize), 0., 1.);
  }
}



void Model::run() {
 // cout << "inside Model::run\n"; //EDIT
  Timer rate_timer, wall_timer;

  //cout << "Created two timer objects, now calling openTrajectoryDCD\n"; //EDIT
  openTrajectoryDCD();
  //cout << "openTrajectoryDCD did not crash the program\n"; //EDIT

  rate_timer.start();
  wall_timer.start();
  cout << " >> Beginning simulation..." << endl;

  while(!done) {
   // //cout << "Inside \"while\" loop, about to start commands\n"; //EDIT
    integrate();
    //cout << "integrate called succesfully\n"; //EDIT
    step++;
    //cout << "step incremented\n"; //EDIT

    if(step % rate_trj == 0) {
      writeCoordinatesDCD();
    }

   if (step % rate_rst == 0){
      writeTimes();
      writeCoordinatesPQR();
      }

    if(step % rate_beta == 0) {
      rate_timer.stop();
      lout << "* Step: " << step << " Walltime: ";
      wall_timer.log_current(&lout);
      lout << " StepRate: " << (rate_timer.duration/rate_beta) << " s/step" << endl;
      printRateConstant();
     
      // check to see if all sessions are done
      bool _done = true;
      for(int i=0; i < sessions.size(); i++) {
        if(! sessions[i]->done) { _done = false; break; }
      }
      done = _done;
      // start timer
      rate_timer.start();
    }
  }

  // Final information
  printRateConstant();
  closeTrajectoryDCD();
}


void Model::printRateConstant() {
  for(int i=0; i < sessions.size(); i++) {
    sessions[i]->printRateConstant();
  }
}


void Model::populateLigands() {
  for(int i=0; i < sessions.size(); i++) {
    sessions[i]->populateLigands();
  }
}


