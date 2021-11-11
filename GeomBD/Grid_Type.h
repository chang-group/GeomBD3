#ifndef _Grid_Type_h_
#define _Grid_Type_h_

#include "Main.h"
#include "Strings.h"
#include "Grid.h"

 
class Grid_Type : public Grid {
  public:
    Grid_Type(string bpm_filename, string atomtype) : Grid(bpm_filename, atomtype) {
    }

};

#endif
