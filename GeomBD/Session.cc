#include "Session.h"
#include "Model.h"
#include "Body.h"
#include "BindingCriteria.h"
#include "Grid_EX.h"
#include "Grid.h"
#include "Strings.h"
#include "Grid_ES.h"
#include "Grid_D.h"
#include "Timer.h"


///Session definition 
Session::Session(Model *m, SimulationConfig s) {
  id = 0;
  model = m;
  type = s;
  Nreplicates = 0;
  Nbind.set_value(0);
  Nexit.set_value(0);
  Ntlim.set_value(0);
  Davg = 0.;
  done = false;
}

Session::~Session() {
}


void Session::populateLigands() {
  for(int i=0; i < Nreplicates; i++) {
    int rndConf = floor(random(0.0, (double)conformations.size())); //Randomly select a ligand conformation
    Body *rc = conformations[rndConf];
    cout << "about to create a new Body\n for replicate " << i+1<< endl;
    Body *bi = new Body(model, this); // Construct new Body object by passing parameters "model" and "this"
    for(int j=0; j < rc->beads.size(); j++) { //Copy over every bead in conformations to bi
      Bead *nb = new Bead();
      Bead *rb = rc->beads[j]; //*rb is pointer to jth bead in the conformation
      nb->R.x = rb->R.x;       //Next, coords and parameters of *rb are copied to *nb
      nb->R.y = rb->R.y;
      nb->R.z = rb->R.z;
      nb->q = rb->q;
      nb->r = rb->r;
      nb->m = rb->m;
      nb->S = rb->S; //New 3/19  
      //cout << "bead " << j+1 << "S assigned to nb: " << nb->S << endl;
      nb->type = rc->beads[j]->type;
      bi->beads.push_back(nb); //*bi is the Body (i.e. ligand) made of all these beads (i.e. atoms). Each Bead nb is pushed into
    }                          //Body bi's bead vector
    bi->define();
    if(model->useRestart) { 
      positionLigandRestart(bi, i); 
      cout << "using rst time from vrt[" << i << "]\n";
      bi->t = model->vrt[i];
      cout<< "t for ligand "<<i+1<< " is " << model->vrt[i] <<endl;
     }
    else { positionLigand(bi); }   
    model->ligands.push_back(bi);   //push ligand bi into vector holding all ligands accross all sessions
    ligands.push_back(bi);          //push ligand into vector holding ligands in one session instance
  }

  Davg = 0.;
  for(int i=0; i < conformations.size(); i++) {
    Davg += conformations[i]->D;
  }
  Davg /= conformations.size();
  cout << " >> Done positioning ligands" << endl;
}

void Session::recordBeta() {
  double _Nbind = (double)Nbind.get_value();
  double _Ndone = (double)Nexit.get_value() + (double)Ntlim.get_value();
  double beta = _Nbind / (_Nbind + _Ndone);
  beta_history.push_back(beta);
  while(beta_history.size() > model->convergence_window) {
    beta_history.pop_front();
  }
}

void Session::recordBindTime(Body *bi) {  //New
  double bt = bi->t;
  bindtimes.push_back(bt);
  model->lout << " ** Pushed back bind time " << bt << endl;
}

double Session::checkConvergence() {
  double m = 0., s = 0., SEM = 0; // New 4/13
  if(model->isNAM or model->logDirect) {  //New if condition, erase if broken
    if(beta_history.size() < model->convergence_window) return -1;

    for(int i=0; i < beta_history.size(); i++) {
      m += beta_history[i];
      }
    m /= beta_history.size();
    if(m == 0) return -1;

    for(int i=0; i < beta_history.size(); i++) {
      s += pow(beta_history[i]-m, 2);
      }
    s = sqrt(s / beta_history.size());
    double soverm = s / m;
    if(soverm <= model->convergence and soverm != 0. and done == false) {
      model->lout << "* Convergence criterion reached. Exiting successfully." << endl;
      done = true;
      for(int i=0; i < ligands.size(); i++) ligands[i]->done = true;
      }

    return soverm;
    }
// Convg of SEM of bind times
  else {
    if(bindtimes.size() < model->convergence_window) return -1;
    for(int i=0; i < bindtimes.size(); i++) {
      m += bindtimes[i];
      }
    m /= bindtimes.size();
    if(m == 0) return -1;

    for(int i=0; i < bindtimes.size(); i++) {
      s += pow(bindtimes[i]-m, 2);
      }
    s = sqrt(s / bindtimes.size());
    SEM = s / sqrt(bindtimes.size());   
    double SEMoverm = SEM / m;
    if(SEMoverm <= model->/*SEM_*/convergence and SEMoverm != 0. and done == false) {
      model->lout << "* Convergence criterion reached. Exiting successfully." << endl;
      done = true;
      for(int i=0; i < ligands.size(); i++) ligands[i]->done = true;
      }
    return SEMoverm;
    }
// END new convg...
}

//SessionA definition
SessionA::SessionA(Model *m) : Session(m, CONFIGURATION_RADIAL) {
  b = 0.;
  q = 0.;
  q2 = 0.;
}

void SessionA::positionLigandRestart(Body *bi, int iter) {
 vertex Q;
 int bcnt=0;
 double tqx=0., tqy=0., tqz=0., tm=0.;
  bi->center(); // Center before changing R.x, etc., center() applies an equal/opp vector to Body's COM
  for (int r=bi->beads.size()*3*iter; r<bi->beads.size()*3*(iter+1); r+=3){
     bi->beads[bcnt]->R.x = model->vrc[r];
     bi->beads[bcnt]->R.y = model->vrc[r+1];
     bi->beads[bcnt]->R.z = model->vrc[r+2];
     tqx += model->vrc[r]  *  bi->beads[bcnt]->m ;
     tqy += model->vrc[r+1] * bi->beads[bcnt]->m ;
     tqz += model->vrc[r+2] * bi->beads[bcnt]->m ;
     tm += bi->beads[bcnt]->m;
     bcnt++;
    }
  Q.x = tqx/tm; Q.y=tqy/tm; Q.z=tqz/tm;
  cout<<"ligdn COM is " <<Q.x<<" "<<Q.y<<" "<<Q.z<<endl;
  bi->R.x = Q.x;
  bi->R.y = Q.y; 
  bi->R.z = Q.z;
  cout << "just positioned a ligand inside positionLigandRestart()\n";
}

void SessionA::positionLigand(Body *bi) {
 if (model->rand_start == 0 and b > model->wall and model->usePlane){ 
  cout << " ! Error: Ligand starting positions are outside of simulation boundary! (must have b <= exbound)" << endl; 
  std::exit(0); 
  }

vertex Q = { random(-1., 1.), random(-1., 1.), random(-1., 1.) };
// Start ligands from single point
 if (model->point_start) {
  Q.x = model->pnt_start.x;
  Q.y = model->pnt_start.y;
  Q.z = model->pnt_start.z;
 }
//Square plane
 else {
  double pw = model->planeWidth, wall = model->wall;
  if (model->boundary == 1) {
   Q.x *= pw;
   Q.y *= pw;
   Q.z = model->bPlnHgt;

   if (model->rand_start == 1) {
     bool penetrating = true;
      while (penetrating) {
       Q.x  = random(-1*wall, wall);
       Q.y  = random(-1*wall, wall);
       Q.z  = random(0., model->ceiling); // Change later !!! 
       if (!bi->checkPenetration(Q.x, Q.y, Q.z)) penetrating=false;
      }
    }
  }
//Hexagonal plane
 else if (model->boundary == 2) {
  Q.y *= pw; 
  if (Q.y >= 0) {Q.x *= (pw-Q.y)*tan(M_PI/6) + pw*tan(M_PI/6);}
      else      {Q.x *= (-pw-Q.y)*tan(M_PI/6) - pw*tan(M_PI/6);}
  Q.z = model->bPlnHgt;

    if (model->rand_start == 1) {
      bool penetrating = true;
      while (penetrating) {
      Q.y = random(-1*wall, wall);
      if (Q.y >= 0){ Q.x = random(-1 * ((wall-Q.y)*tan(M_PI/6) + wall*tan(M_PI/6)), ((wall-Q.y)*tan(M_PI/6) + wall*tan(M_PI/6))); }
      else        { Q.x = random(-1 * ((-wall - Q.y)*tan(M_PI/6) - wall*tan(M_PI/6)), ((-wall - Q.y)*tan(M_PI/6) - wall*tan(M_PI/6))); }
      Q.z = random(0., model->ceiling);
      if (!bi->checkPenetration(Q.x, Q.y, Q.z)) penetrating=false; 
      }
    }
  }
//NAM b sphere
 else if (model->boundary == 3 /*or b != 0.*/) { //New b start for dir/indir
  double l = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
  Q.x *= b / l;
  Q.y *= b / l;
  Q.z *= b / l;
  }
//Spherical boundary
 else if (model->boundary == 4) {
   if (model->rand_start == 1) {
    bool penetrating = true;
    while (penetrating) {
      bool out_of_bounds=true;
      while(out_of_bounds){
        Q.x = random(-wall, wall); 
        Q.y = random(-wall, wall);
        Q.z = random(-wall, wall);
        if( (Q.x*Q.x + Q.y*Q.y + Q.z*Q.z) < wall*wall ) out_of_bounds=false;
        }
      if (!bi->checkPenetration(Q.x, Q.y, Q.z)) penetrating=false; 
      }
    }
  }
 }

  if(!model->rand_start==1) {
    bi->center();
    if (!model->point_start) { 
      if(model->boundary==1 or model->boundary==2) bi->translate(model->center.x + Q.x, model->center.y + Q.y, model->receptor_min.z + Q.z);
      if(model->boundary==3 or model->boundary==4) bi->translate(model->center.x + Q.x, model->center.y + Q.y, model->center.z + Q.z);
      bi->rotate(random(0., M_PI), random(0., M_PI), random(0., M_PI)); 
      }
    else bi->translate(Q.x, Q.y, Q.z);
    }   
}

void SessionA::positionIndirect(Body *bi) { //Called upon first crossing of dircut in dir/indir sim
  vertex Q = { random(-1., 1.), random(-1., 1.), random(-1., 1.) };
  double l = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
  Q.x *= b / l;
  Q.y *= b / l;
  Q.z *= b / l;
  bi->center();
  bi->rotate(random(0., M_PI), random(0., M_PI), random(0., M_PI));
  bi->translate(model->center.x + Q.x, model->center.y + Q.y, model->center.z + Q.z);
}

void SessionA::printRateConstant() {
  int Ndone = Nbind.get_value() + Nexit.get_value();
  double conv = checkConvergence(); 
  if(Ndone == 0) {
    double kb = 4. * M_PI * b * Davg   *Na*1e12*1e-27;
    model->lout << "   (session " << id << ") ";
    model->lout << "kd(b)=" << kb << " ";
    model->lout << "b=" << b << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
    return;
  }

  for(int bsi=0; bsi < bindingCriteria.size(); bsi++) {
    double tavg = bindingCriteria[bsi]->t_avgt.get_value() / bindingCriteria[bsi]->Nbind.get_value();
    double Ndir = bindingCriteria[bsi]->dirTimes.size(), Nindir = bindingCriteria[bsi]->indirTimes.size();
    double k, kb, B;
    if(model->logDirect) B = Nindir / (Nindir + double(Nexit.get_value())); //Denominator doesn't include direct binders; they are separate
    else B = ((double)bindingCriteria[bsi]->Nbind.get_value()) / ((double)Ndone); 
    if(model->isNAM) kb = 4. * M_PI * b * Davg;
    else kb = 4. * M_PI * model->dircut * Davg; //New 4/28
    if(model->isNAM or model->logDirect) k = (kb * B) / (1 - ((1 - B)*b/q));
    else k = model->volume / tavg; //New
    k *= Na * 1e12 * 1e-27;
    model->lout << "   (session " << id << " bs " << bsi << ") ";
    model->lout << "k_on = " << k << " M⁻¹s⁻¹ ";
    if(model->logDirect) { // New 4/16
       double Bdir = Ndir / (Ndir + Nindir + (double)Nexit.get_value());
       double sum_t=0.;
       for (auto& ts : bindingCriteria[bsi]->dirTimes) sum_t += ts; //Sums over all elements in vector
       double tavg_dir = sum_t / Ndir;
       model->lout << "tavg_dir=" << tavg_dir << " Bdir=" << Bdir << ", k_dir = " << Bdir / (tavg_dir * 1e-12) << " " << endl;
       model->lout << "Ndirect=" << bindingCriteria[bsi]->dirTimes.size() << " ";
       model->lout << "Nindirect=" << bindingCriteria[bsi]->indirTimes.size() << " ";
       }
    model->lout << "Nbind=" << bindingCriteria[bsi]->Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    if(model->isNAM or model->logDirect) model->lout << "β=" << B << " conv=" << conv << " ";
    else model->lout << "SEM=" << SEM << " ps conv=" << conv << " ";
    if(done) model->lout << "(convergence reached) ";
    model->lout << "kd(b)=" << kb * Na * 1e12 * 1e-27 << " ";
    model->lout << "b=" << b << " " << "q=" << q << " " << "Davg=" << Davg << " ";
    //model->lout << "Boundary type = " << model->boundary << endl; //New
  }
 
  double B = ((double)Nbind.get_value()) / ((double)Ndone);
  if(bindingCriteria.size() > 1) {
    double tavg = t_avgt.get_value() / Nbind.get_value();
    double kb = 4. * M_PI * b * Davg;
    double k = (kb * B) / (1 - ((1 - B)*b/q));  
    //double k = 50000000. / tavg;    //New k with no exit criteria, beta not needed
    k *= Na * 1e12 * 1e-27;
    model->lout << "   (session " << id << ") ";
    model->lout << "k_on = " << k << " M⁻¹s⁻¹ ";  
    model->lout << "Nbind=" << Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "β=" << B << " ";
    model->lout << "kd(b)=" << kb * Na * 1e12 * 1e-27 << " ";
    model->lout << "b=" << b << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
  }

  if(model->max_simulations > 0 and Ndone >= model->max_simulations) {
    model->lout << "> Maximum simulations reached. Exiting." << endl;
    model->done = true;
  }
}


void SessionA::checkLigand(Body *bi) {
// Only check q sphere penetration for spherical NAM or direct/indirect transfer simulation
 if (model->boundary == 3){
  if(model->logDirect) q = model->wall; //New 4/16 set q here since not set in Model_Input
  double dr[3], l2;
  dr[0] = bi->R.x - model->center.x;
  dr[1] = bi->R.y - model->center.y;
  dr[2] = bi->R.z - model->center.z;
  l2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
  if(l2 >= q2) {
    if(model->logExiters)
      model->lout << "#" << id << "\t Escape event at t=" << bi->t << " ps  (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->exited = true;
    bi->session->positionLigand(bi);
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
    bi->indirect = false; //New 4/15
  }
 }
}
