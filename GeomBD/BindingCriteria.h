#ifndef BindingCriteria_h_
#define BindingCriteria_h_
#include <cilk/reducer_opadd.h>
#include <iostream>

struct BindingPair { //The basic unit of a binding criterion: a distance, an xyz coord, and a lig atom #
  double r2;
  vertex target;
  int laid;
};

class BindingCriteria {
  public:
    vector< vector<BindingPair> > pairs; //Holds all the bindingSets in simulation
    vector<double> critNbind; 
    vector<double> critTimes;
    vector<double> dirTimes;
    vector<double> indirTimes;
    cilk::reducer< cilk::op_add<int> > Nbind;
    cilk::reducer< cilk::op_add<double> > t_avgt;
    bool AND;
    bool ANDOR; //EDIT 6/16/18 @ 12:32pm
    int criterionID; // New 3/8/21


  public:
    BindingCriteria(bool cAND=true, bool cANDOR=false) {  //"bindgroups" in input file makes AND false, ANDOR true
      Nbind.set_value(0);
      t_avgt.set_value(0.);
      AND = cAND;
      ANDOR = cANDOR;
      vector<BindingPair> bindingSet; //Each bindingSet holds 1 or more binding pairs (> 1 pair for AND conditions)
      pairs.push_back(bindingSet);
    }

    ~BindingCriteria() {}

    //EDIT- Modify addpair by adding in parameter "grpID" determining the AND group that the BindingPair will be added to
    void addPair(double rx, double ry, double rz, int lig_id, double r, int grpID=0) {
      if(pairs.size() == grpID){ 
        vector<BindingPair> bindingSet;
        pairs.push_back(bindingSet);
      }
      BindingPair bp;
      bp.r2 = r*r;
      bp.target.x = rx;
      bp.target.y = ry;
      bp.target.z = rz;
      bp.laid = lig_id;
      pairs[grpID].push_back(bp); //push bp into bindignSet[grpID]
    }
    //ANDOR if condition checks for any true AND grouping, if found returns true for binding.
    bool checkBinding(Body *ligand) {
      if(ANDOR){
        for(int grpID=0; grpID < pairs.size(); grpID++){ // Loops over all vectors in "pairs"
          if(checkBindingAND(ligand, grpID))
            return true;
        }
      }
      else if(AND) return checkBindingAND(ligand);
      return checkBindingOR(ligand); //checkBindingOR only called if checkBindingAND does not return true above
    }
    
    bool checkBindingAND(Body *ligand, int grpID=0) {
      for(int i=0; i < pairs[grpID].size(); i++) { //Loops over all criteria in this bindingSet
        vertex dr;
        if(pairs[grpID][i].laid < 0) { //This seems to do a check based on ligand COM if laid is -1 I guess?
          dr.x = ligand->R.x - pairs[grpID][i].target.x;
          dr.y = ligand->R.y - pairs[grpID][i].target.y;
          dr.z = ligand->R.z - pairs[grpID][i].target.z;
        } else {
          Bead *lb = ligand->beads[pairs[grpID][i].laid];
          dr.x = lb->R.x - pairs[grpID][i].target.x;
          dr.y = lb->R.y - pairs[grpID][i].target.y;
          dr.z = lb->R.z - pairs[grpID][i].target.z;
        }
        double l2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
        if(l2 > pairs[grpID][i].r2) return false;  //Since this is AND grouping, as soon as find one pair that is not bound, return FALSE
      }
      criterionID = grpID; //New 3/8/21
      critNbind[grpID]++; // New
      critTimes[grpID] += ligand->t; //New
      return true;  //Otherwise, return TRUE because we went thru for loop w/o finding any that weren't bound TC
    }

    bool checkBindingOR(Body *ligand) {
        for(int i=0; i < pairs.size(); i++){
          vertex dr;
          if(pairs[i][0].laid < 0) {
            dr.x = ligand->R.x - pairs[i][0].target.x;
            dr.y = ligand->R.y - pairs[i][0].target.y;
            dr.z = ligand->R.z - pairs[i][0].target.z;
          } else {
            Bead *lb = ligand->beads[pairs[i][0].laid];
            dr.x = lb->R.x - pairs[i][0].target.x;
            dr.y = lb->R.y - pairs[i][0].target.y;
            dr.z = lb->R.z - pairs[i][0].target.z;
          }
          double l2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
          if(l2 <= pairs[i][0].r2) {
             criterionID = i; //New 3/8/21
             critNbind[i]++; // New
             critTimes[i] += ligand->t; //New
             return true;  // if any OR conditiion is met, return TRUE and stop checking TC
            }
        }
      return false;
    }


};


#endif
