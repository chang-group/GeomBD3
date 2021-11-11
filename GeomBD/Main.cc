
#include "Main.h"
#include "Timer.h"
#include "Model.h"


bool getInputWithFlag(int argc, char **argv, char flag, string *value) 
{
  int  opt;                                                                                                                                                                                                     
  bool bopt = false;
  char gopt[2] = { flag, ':' };

  for(int i=1; i < argc; i++) 
  {
    if(argv[i][0] == '-' && argv[i][1] == flag) 
    {
      if(argv[i+1][0] != '-') 
      {
        *value = argv[i+1];
        bopt = true;
        i++;
        break;
      }
    }
  }

  return bopt;
}


void usage() 
{
  printf("Usage: GeomBD3 -i INPUTFILE -o TRAJECTORY.pqr -l LOGFILE.log (Optional: -c RESTART_COOR.crd -t RESTART_TIMERS.t)\n");
}


static Model *model = NULL;

void term(int signal) 
{
  if(model) 
  {
    model->done = true;
  }
}

int main(int argc, char **argv) 
{
  string stoken; 
  string fldfn;       //input file name
  string trjfn;       //trajectory file name
  string logfn;       //log file name
  string rst_t, rst_crd;  //New
  srand(time(NULL));

  if(!getInputWithFlag(argc, argv, 'i', &fldfn)) 
    { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'o', &trjfn)) 
    { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'l', &logfn)) 
    { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'c', &rst_crd)) //New
    { cout << "No restart coordinates given\n"; }
  if(!getInputWithFlag(argc, argv, 't', &rst_t)) //New
    { cout << "No restart timers given\n"; }

  // Exit gracefully if possible
  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = term;
  sigaction(SIGINT, &action, NULL);
  sigaction(SIGTERM, &action, NULL);
  sigaction(SIGQUIT, &action, NULL);

  // Create model
  //if(getInputWithFlag(argc, argv, 'c', &rst_crd)) {
    // code for new model constructor w 5 agrs
    model = new Model(fldfn, trjfn, logfn, rst_crd, rst_t);
  //}
  //else { model = new Model(fldfn, trjfn, logfn); } // These file names become ifn, ofn, and lfn by the Model constructor

  model->lout << "I CAN use classes too, I swear!\n";
  model->lout << "* Running on " << model->threads << " threads" << endl;
  Timer *timer = new Timer();
  model->lout << "We created a timer\n";
  timer->start();
  model->lout << "Timer started succesfully\nNow calling Model::run()\n";
  model->run();
  //cout << " Called model->run()\n";
  // Above and below this line is where simulation starts then stops
  model->lout << "Model started running, now calling Model::stop()\n";
  timer->stop();
  model->lout << "The model managed to stop\n";
  timer->print(&model->lout);


  delete timer;
  delete model;

  return EXIT_SUCCESS;
}



