#ifndef _Grid_D_h_
#define _Grid_D_h_

#include "Main.h"
#include "Strings.h"
#include "Grid.h"

 
class Grid_D : public Grid {
  public:
    Grid_D(string bpm_filename, string atomtype) : Grid(bpm_filename, atomtype) {
    }


};

#endif
