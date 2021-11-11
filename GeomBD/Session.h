#ifndef Session_h_
#define Session_h_

#include "Main.h"

class Model;
class Body;
class BindingCriteria;

 
class Session {
  friend class Model;

  protected:
    int id;
    Model *model;
    SimulationConfig type;
    int Nreplicates;
    cilk::reducer< cilk::op_add<int> > Nbind;
    cilk::reducer< cilk::op_add<int> > Nexit;
    cilk::reducer< cilk::op_add<int> > Ntlim;
    cilk::reducer< cilk::op_add<double> > t_avgt;
    double Davg;
    vector< Body* > conformations;
    vector< Body* > ligands;
    vector< BindingCriteria* > bindingCriteria;
    deque<double> beta_history;
    vector<double> bindtimes; // New 3/4/21
    bool done;

  public:
    Session(Model *m, SimulationConfig s);
    ~Session();
     
    virtual void populateLigands();
    virtual void positionLigand(Body *body) { };
    virtual void positionLigandRestart(Body *body, int it) { }; // New
    virtual void positionIndirect(Body *body) { }; //New 5/18
    virtual void printRateConstant() { };
    virtual void checkLigand(Body *body) { };
    virtual void recordBeta();
    virtual void recordBindTime(Body *body); //New 3/4/21
    virtual double checkConvergence();
    virtual void finalize() { };
    double SEM;
};

//All simulations happen in one session type called SessionA 05/27/21
class SessionA : public Session {
  friend class Model;
  protected:
    double b;
    double q, q2;

  public:
    SessionA(Model *m);

    virtual void positionLigand(Body *body);
    virtual void positionLigandRestart(Body *body, int it); //New 12/4
    virtual void positionIndirect(Body *body); //New 5/18
    virtual void printRateConstant();
    virtual void checkLigand(Body *body);

};

//SessionDirect depricated 05/27/21

/*
class SessionDirect : public Session {
  friend class Model;
  protected:
    vertex start;
    double b;
    double q, q2;

  public:
    SessionDirect(Model *m);

    virtual void positionLigand(Body *body);
    virtual void printRateConstant();
    virtual void checkLigand(Body *body);

};
*/

#endif
