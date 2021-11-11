#ifndef Main_h_
#define Main_h_

#include <cstdlib>
#include <cstdio>
#include <csignal>
#include <memory.h>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <algorithm> 
#include <sys/time.h>
#include <sys/stat.h>
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
#include "cilk/reducer_opadd.h"
#include "cilk/reducer_list.h"
#include <mkl.h>
#include <limits.h>

using namespace std;


/*
 * Constants
 */
#define kB 1.9858775e-3 // kcal/molK
#define kC 332.0623183 // coulomb constant
#define Na 6.0221415e23
#define LperA3 1.0e-27
#define nMperM 1.0e9


/*
 * Conversions
 */
#define CONCENTRATION(Nparticles, BOXL) (((()Nparticles / Na) / ((BOXL * BOXL * BOXL) * LperA3)) * nMperM)
#define BOXSIZEFORnMCONC(nMCONC) pow((nMperM * (1.0 / Na) / (nMCONC * LperA3)), 1./3.)

/*
 * Convenience
 */

inline double random(double rangeStart, double rangeEnd) {
  return (((double)rand() / (double)INT_MAX) * (rangeEnd - rangeStart)) + rangeStart;
}

inline bool file_exists(const string& fn) {
  struct stat buffer;   
  return (stat (fn.c_str(), &buffer) == 0); 
}


/*
 * Type definitions
 */

struct vertex {
  double x;
  double y;
  double z;
};

inline double vertex_sqmagnitude(vertex &v) {
  return (v.x*v.x + v.y*v.y + v.z*v.z);
}

inline double vertex_magnitude(vertex &v) {
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

enum SimulationConfig {
  CONFIGURATION_RADIAL,
  CONFIGURATION_ABSOLUTE_RADIAL,
  CONFIGURATION_ABSOLUTE_PERIODIC,
  CONFIGURATION_MILESTONE
};

/*
struct BindingSite {
  double r2;
  double x, y, z;
};
*/


typedef unsigned int unt;

#endif
