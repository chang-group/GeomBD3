#ifndef _Grid_ES_h_
#define _Grid_ES_h_

#include "Main.h"
#include "Strings.h"
#include "Grid.h"

 
class Grid_ES : public Grid {
  public:
    Grid_ES(string bpm_filename, string atomtype) : Grid(bpm_filename, atomtype) {
    }

};

#endif
