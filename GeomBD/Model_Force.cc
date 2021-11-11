
#include "Model.h"
#include "Body.h"
#include "BindingCriteria.h"



void Model::integrate() {

  generateNormal();

  cilk_for(int il=0; il < ligands.size(); il++) {
    double oe_elec, oe_d, oe_lj, toe_elec, toe_d, toe_lj; //New output E variables
    Body *Bi = ligands[il];
    if(! Bi->done) {
      bool onGrid = false, associated = false;
      double E;

      // Calculate forces
      for(int i=0; i < Bi->beads.size(); i++) {
        Bead *bi = Bi->beads[i];
        bi->F.x = 0.;
        bi->F.y = 0.;
        bi->F.z = 0.;
        vertex dF;
        dF.x = 0.;
        dF.y = 0.;
        dF.z = 0.;
        
        if(bi->q != 0.) {
          for(int es=0; es < esmaps.size(); es++) {
            if(esmaps[es]->approximate_force(&bi->R, &dF, &E, fd_order, bi->q)) {
              onGrid = true;
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
              if(fabs(E) > 0 and il == 0 and step > 0 and step % rate_trj == 0) { 
                esmaps[es]->approximate_potential(&bi->R, &oe_elec); //New
                oe_elec *= bi->q;
                toe_elec += oe_elec;
              }
            }
          }
        }
          for(int ds=0; ds < dmaps.size(); ds++) {
            if(dmaps[ds]->approximate_force(&bi->R, &dF, &E, fd_order, bi->S)) { //New 3/22
              onGrid = true;
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
              if(fabs(E) > 0) {
                associated = true;
                if(il == 0 and step > 0 and step % rate_trj == 0) {
                  dmaps[ds]->approximate_potential(&bi->R, &oe_d); //New
                  oe_d *= bi->S;
                  toe_d += oe_d;
                }
              }
            }
          }
        

        for(int tmap=0; tmap < typemaps.size(); tmap++) {
          if(typemaps[tmap]->type == bi->type) { //Only check map is bead type matches map type
            if(typemaps[tmap]->approximate_force(&bi->R, &dF, &E, fd_order, 1.0)) {
              onGrid = true;
              if(fabs(E) > 0) {
                associated = true;
                if(il == 0 and step > 0 and step % rate_trj == 0) {
                  typemaps[tmap]->approximate_potential(&bi->R, &oe_lj); //New
                  //oe_lj *= bi->S;
                  toe_lj += oe_lj;
                }
              }
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
            }
          }
        }
      } //End loop over atoms

      if(step>0 and step%rate_trj == 0 and il==0) lout << " >> Replicate 1 U_elec = " << toe_elec << " kcal/mol, U_d = " << toe_d << "kcal/mol, U_lj = " << toe_lj << " kcal/mol" << endl;

      // Calculate force magnitude
      vertex Rcom = { Bi->R.x - center.x, Bi->R.y - center.y, Bi->R.z - center.z }; //Is this distance of lig COM from rec center? I think
      double Rcom_mag = sqrt(Rcom.x*Rcom.x + Rcom.y*Rcom.y + Rcom.z*Rcom.z);
      if(Rcom_mag != 0.) {
        Rcom.x /= Rcom_mag;
        Rcom.y /= Rcom_mag;
        Rcom.z /= Rcom_mag;
      }
      vertex Fi = Bi->F; //has Bi F been changed before this point? Or maybe its whatever it was in prev step, Fi may be force initial
      double Fi_mag = sqrt(Fi.x*Fi.x + Fi.y*Fi.y + Fi.z*Fi.z);
      if(Fi_mag != 0.) {
        Fi.x /= Fi_mag;
        Fi.y /= Fi_mag;
        Fi.z /= Fi_mag;
      }
      double Rcom_dot_Fi = Rcom.x*Fi.x + Rcom.y*Fi.y + Rcom.z*Fi.z;
      double New_mF = Rcom_dot_Fi * Fi_mag; //What does this do?

      // Propogate bead forces to body
      Bi->F.x = 0.;
      Bi->F.y = 0.;
      Bi->F.z = 0.;
      Bi->Fa.x = 0.;
      Bi->Fa.y = 0.;
      Bi->Fa.z = 0.;

      for(int k=0; k < Bi->beads.size(); k++) {
        Bead *bk = Bi->beads[k];

        if(bk->F.x == 0. and bk->F.y == 0. and bk->F.z == 0.) continue;

        Bi->F.x += bk->F.x;      //adds up force from all atoms in body to get Bi->F.x, total x force, so on for y z
        Bi->F.y += bk->F.y;
        Bi->F.z += bk->F.z;

        vertex A = { bk->R.x - Bi->R.x, bk->R.y - Bi->R.y, bk->R.z - Bi->R.z }, B; //Finding bead k's distance from body COM to calc rot Forces
        B.x = A.y * bk->F.z - A.z * bk->F.y;     // y axis * z Force - z axis * y force give rot about x axis, etc
        B.y = A.z * bk->F.x - A.x * bk->F.z;    
        B.z = A.x * bk->F.y - A.y * bk->F.x;

        Bi->Fa.x += B.x;
        Bi->Fa.y += B.y;
        Bi->Fa.z += B.z;
      }

      // Determine timestep
      double dr[3], radius2, radius;
      dr[0] = Bi->R.x - center.x;
      dr[1] = Bi->R.y - center.y;
      dr[2] = Bi->R.z - center.z;
      radius2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

      if(radius2 > dt_scale_start and radius2 < dt_scale_end) {
        double s = (radius2 - dt_scale_start) / (dt_scale_end - dt_scale_start);
        Bi->dt = dt_fine + s * (dt_coarse - dt_fine);
      } else {
        if(radius2 <= dt_scale_start) {
          Bi->dt = dt_fine;
          if(New_mF < -1.0) {
            Bi->dt /= (New_mF / -1.0);
          }
        }
        if(radius2 >= dt_scale_end) Bi->dt = dt_coarse;
      }

      // Set direct/indirect label
      if(logDirect) {
        if(!Bi->indirect){ //Only set label and timer at first crossing
           if(sqrt(radius2) > dircut) { Bi->indirect=true; Bi->t = 0.; Bi->session->positionIndirect(Bi); } //New 05/18
           }
        //if(step % rate_trj ==0 and il==0) lout << " ** radius = "<< sqrt(radius2) << ", cut = " << dircut <<", dir/indir = " << Bi->indirect << endl;
        }
      // Backup coordinates in case of interpenetration
      if(exmaps.size() > 0)
        Bi->save();

      // Integrate
      vertex *Si = &rand[il];
      vertex *Sai = &rand[ligands.size()+il];
 
      if(step>0 and step%rate_trj == 0 and il==0) lout << " * F xyz on ligand 1: " << Bi->F.x <<" "<<Bi->F.y<<" "<<Bi->F.z<<"\n * Random F xyz: " <<Si->x<<" "<<Si->y<<" "<<Si->z<<endl;
     
      vertex dR, dRa;
      double A = sqrt(2. * Bi->D * Bi->dt); 
      double dtOVERkBT = Bi->dt / (kB * T);
      double B = Bi->D * dtOVERkBT;
      if(flux != 0) dR.x = (A * Si->x) + (B * (Bi->F.x + flux) ); //New
      else dR.x = (A * Si->x) + (B * Bi->F.x);
      dR.y = (A * Si->y) + (B * Bi->F.y); 
      dR.z = (A * Si->z) + (B * Bi->F.z);
      Bi->translate(dR.x, dR.y, dR.z);

      double C = sqrt(2 * Bi->Da * Bi->dt);
      double D = Bi->Da * dtOVERkBT;
      dRa.x = (C * Sai->x) + (D * Bi->Fa.x);
      dRa.y = (C * Sai->y) + (D * Bi->Fa.y);
      dRa.z = (C * Sai->z) + (D * Bi->Fa.z);
      Bi->rotate(dRa.x, dRa.y, dRa.z);

     //Cuboid periodic translation
     if (boundary == 1){
      double lx = center.x + wall;
      double ly = center.y + wall;
      double x;
      double y;
      x = (Bi->R.x);
      y = (Bi->R.y);
      if (x > lx)          { Bi->translate(-2.*wall, 0., 0.); }
      if (x < lx-(2*wall)) { Bi->translate(2.*wall, 0., 0.);  }   
      if (y > ly)          { Bi->translate(0., -2.*wall, 0.); }
      if (y < ly-(2*wall)) { Bi->translate(0., 2.*wall, 0.);  }
      if (Bi->R.z < receptor_min.z) { Bi->restore();}
      if (Bi->R.z > receptor_min.z + ceiling) { Bi->restore();}
      }
     //Hexagonal periodic translation
     if (boundary == 2){
      double hx1  =  center.x + wall*(tan(M_PI/6)) + ( (wall-(Bi->R.y-center.y))*(tan(M_PI/6)) );
      double hx2  =  center.x - wall*(tan(M_PI/6)) - ( (wall-(Bi->R.y-center.y))*(tan(M_PI/6)) );
      double nhx1 =  center.x + wall*(tan(M_PI/6)) + ( (wall-(center.y-Bi->R.y))*(tan(M_PI/6)) ); 
      double nhx2 =  center.x - wall*(tan(M_PI/6)) - ( (wall-(center.y-Bi->R.y))*(tan(M_PI/6)) );
      double hy = center.y + wall;              
      double x, y;
      x = (Bi->R.x);
      y = (Bi->R.y);
      if(x >= hx1 && y >= center.y)  { Bi->translate((-1.0)*(sqrt(3.0*wall*wall)), -wall, 0.); }
      if(x <  hx2 && y >= center.y)  { Bi->translate(sqrt(3.0*wall*wall), -wall, 0.); }
      if(x >= nhx1 && y <= center.y) { Bi->translate((-1.0)*(sqrt(3.0*wall*wall)), wall, 0.); }
      if(x <  nhx2 && y <= center.y) { Bi->translate(sqrt(3.0*wall*wall), wall, 0.); }
      if(y >= hy)                    { Bi->translate(0., -2.0*wall, 0.); }
      if(y < hy-(2*wall))            { Bi->translate(0., 2.0*wall, 0.); }
      if(Bi->R.z < receptor_min.z)   { Bi->restore();}
      if(Bi->R.z > receptor_min.z + ceiling) { Bi->restore();}
      }
      //Hard spherical boundary
      if (boundary == 4){
       double l = pow(Bi->R.x-center.x, 2) + pow(Bi->R.y-center.y, 2) + pow(Bi->R.z-center.z, 2);
       if (l >= wall*wall) Bi->restore();
       }

      // Exclusion check
      bool penetrating = false;
              
      for(int i=0; i < Bi->beads.size(); i++) {
        Bead *bi = Bi->beads[i];
        for(int ex=0; ex < exmaps.size(); ex++) {
          if(exmaps[ex]->value(&bi->R) > 0) {
            penetrating = true;
          }
          if(penetrating) Bi->restore();
            break;
         }
        } 

     // Debug check
      if(debug_map) {
        bool debug_penetrate = false;
        for(int i=0; i < Bi->beads.size(); i++) {
          Bead *bi = Bi->beads[i];
          if(debug_map->value(&bi->R) > 0) {
            debug_penetrate = true;
          }
          if(debug_penetrate) break;
        }
        if(debug_penetrate) {
          Bi->mF = New_mF;
        } else {
          Bi->mF = 0.;
        }
      }

      // Increment time, record dwell-time
      Bi->t += Bi->dt;

      if(onGrid and associated) {
        Bi->t_dwell += Bi->dt;
        Bi->t_dwell_total += Bi->dt;
        if(Bi->t_dwell > Bi->t_dwell_max) Bi->t_dwell_max = Bi->t_dwell;
      } else {                                                           //Else handles ligs that assoc then dissoc
        if(Bi->t_dwell > 0.) {                                           
          if(Bi->t_dwell > Bi->t_dwell_max) Bi->t_dwell_max = Bi->t_dwell;
          Bi->t_dwell = 0.; //Resets if ligand dissociates
        }
      }

      // Session-specific checks
      Bi->session->checkLigand(Bi);

      // Binding criteria check
      for(int bs=0; bs < Bi->session->bindingCriteria.size(); bs++) {
        BindingCriteria* bc = Bi->session->bindingCriteria[bs];
        if(bc->checkBinding(Bi)) {
          if(writeBinders) {  // Should we write the bound conformation?
            fstream boutf(bfn.c_str(), ios::out | ios::app);
            Bi->writePDB(boutf, 'A');
            boutf << "TER" << endl;
            boutf.close();
          }
          if(logBinders) {
            lout << "#" << Bi->session->id << "\t Binding criterion " << bc->criterionID+1 << " satisfied at t=" << Bi->t
                 << " ps  (t_dwell=" << Bi->t_dwell << "ps, max=" << Bi->t_dwell_max << "ps, total=" << Bi->t_dwell_total << "ps, \n";  
            lout << "\tBind Times:\n";
            for(int i=0; i<bc->pairs.size(); i++) {
              if(!bc->critNbind[i]==0)
                lout << "\t Criterion " << i+1 <<": Nbind="<<bc->critNbind[i] << " tavg=" << bc->critTimes[i]/bc->critNbind[i] << endl;
              }
            }
          Bi->bound = true;
          *Bi->session->t_avgt += Bi->t;
          *bc->Nbind += 1;
          *bc->t_avgt += Bi->t;
          if(Bi->indirect) bc->indirTimes.push_back(Bi->t); //New
          else bc->dirTimes.push_back(Bi->t);               //New
          Bi->session->recordBindTime(Bi);
          Bi->session->positionLigand(Bi);
          Bi->t = 0.;
          Bi->t_dwell = 0.;
          Bi->t_dwell_max = 0.;
          Bi->t_dwell_total = 0.;
          Bi->indirect = false; //New 4/15
        }
      }

    }
  }

  cilk_sync;

  for(int il=0; il < ligands.size(); il++) {
    Body *Bi = ligands[il];
    if(Bi->bound == true) {
      *Bi->session->Nbind += 1;
      Bi->session->recordBeta();
      Bi->bound = false;
    }
    if(Bi->exited == true) {
      *Bi->session->Nexit += 1;
      Bi->session->recordBeta();
      Bi->exited = false;
    }
    if(Bi->timedout == true) {
      *Bi->session->Nexit += 1;
      Bi->session->recordBeta();
      Bi->timedout = false;
    }
  }
}


