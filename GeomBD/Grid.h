#ifndef _Grid_h_
#define _Grid_h_

#include "Main.h"
#include "Strings.h"
#include "Model.h"
 
class Grid {
  friend class Model;

  protected:
    bool header = true,
         body = false,
         footer = false;

    int N[3];
    int Nt;
    double origin[3];
    double delta;
    double ***data;

  public:
    string type;

  public:
    Grid(string bpm_filename, string atomtype) {
      type = atomtype;

      ifstream fd(bpm_filename.c_str(), ios::in | ios::binary);

      if(!fd.is_open()) {
        cout << "! Failed to load grid file. File not found." << endl;
        exit(-1);
      }

      fd.read((char*)&origin, sizeof(double) * 3);
      fd.read((char*)&N, sizeof(int) * 3);
      fd.read((char*)&delta, sizeof(double));
      Nt = N[0] * N[1] * N[2];

      data = (double***)calloc(N[0], sizeof(double**));
      for(int nx=0; nx < N[0]; nx++) {
        data[nx] = (double**)calloc(N[1], sizeof(double*));              
        for(int ny=0; ny < N[1]; ny++) {
          data[nx][ny] = (double*)calloc(N[2], sizeof(double));              
        }
      }

      for(int nx=0; nx < N[0]; nx++) {
        for(int ny=0; ny < N[1]; ny++) {
          for(int nz=0; nz < N[2]; nz++) {
            fd.read((char*)&data[nx][ny][nz], sizeof(double));
          }
        }
      }

    }

    virtual ~Grid() {
      for(int nx=0; nx < N[0]; nx++) {
        for(int ny=0; ny < N[1]; ny++) {
          free(data[nx][ny]);
        }
        free(data[nx]);
      }
      free(data);
    }

    bool coordinateToGrid(vertex *R, int *G) {
      G[0] = (int)floor(((R->x-origin[0])/delta) + 0.5);
      if(G[0] <= 0 or G[0] >= N[0]-1) {
        return false;
      }
      G[1] = (int)floor(((R->y-origin[1])/delta) + 0.5);
      if(G[1] <= 0 or G[1] >= N[1]-1) {
        return false;
      }
      G[2] = (int)floor(((R->z-origin[2])/delta) + 0.5);
      if(G[2] <= 0 or G[2] >= N[2]-1) {
        return false;
      }

      return true;
    }

    bool coordinateToGrid(vertex *R, int *G, float *Gf) {
      Gf[0] = (R->x - origin[0]) / delta;
      G[0] = (int)floor(Gf[0]);
      if(G[0] <= 0 or G[0] >= N[0]-1) {
        return false;
      }
      Gf[1] = (R->y - origin[1]) / delta;
      G[1] = (int)floor(Gf[1]);
      if(G[1] <= 0 or G[1] >= N[1]-1) {
        return false;
      }
      Gf[2] = (R->z - origin[2]) / delta;
      G[2] = (int)floor(Gf[2]);
      if(G[2] <= 0 or G[2] >= N[2]-1) {
        return false;
      }

      return true;
    }


    /* trilinear interpolation */
    bool approximate_potential(vertex *R, double *e, double scale=1.0) {
      int grid[3];
      float fgrid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid, (float*)&fgrid);
      if(!onGrid) return false;

      double tx = fgrid[0] - grid[0];   //Does this find how far the ligand is from an acutal grid point?
      double ty = fgrid[1] - grid[1];
      double tz = fgrid[2] - grid[2];

      double c000 = data[grid[0]][grid[1]][grid[2]];
      double c100 = data[grid[0]+1][grid[1]][grid[2]];
      double c010 = data[grid[0]][grid[1]+1][grid[2]];
      double c001 = data[grid[0]][grid[1]][grid[2]+1];
      double c110 = data[grid[0]+1][grid[1]+1][grid[2]];
      double c111 = data[grid[0]+1][grid[1]+1][grid[2]+1];
      double c011 = data[grid[0]][grid[1]+1][grid[2]+1];
      double c101 = data[grid[0]+1][grid[1]][grid[2]+1];

      *e = scale * (
           (1.0-tx)*(1.0-ty)*(1.0-tz)*c000 
         + tx*(1.0-ty)*(1.0-tz)*c100 
         + (1.0-tx)*ty*(1.0-tz)*c010 
         + tx*ty*(1.0-tz)*c110 
         + (1.0-tx)*(1.0-ty)*tz*c001 
         + tx*(1.0-ty)*tz*c101 
         + (1.0-tx)*ty*tz*c011
         + tx*ty*tz*c111
      );
      return true;
    }


    bool approximate_force(vertex *R, vertex *F, double *e, int order, double scale) {
      double E, dU[3], te[2], te2[2], te3[2];
      double step = 0.5 * delta;
      double twostep = delta;
      double threestep = 3. * step;
      double fourstep = 4. * step;
      double twelvestep = 12. * step;
      double sixtystep = 60. * step;

      if(order >= 4) {
        cout << "! approximate_force only supports orders 1, 2, and 3" << endl;
        return false;
      }


      /* obtain trilinear intropolation of potential */
      if(! approximate_potential(R, &E)) return false;

      /* Using central difference with first/second/third order approximation to first derivative */
      /* x axis */
      vertex e_v = *R;
      e_v.x += step; //find potential forward one step
      if(! approximate_potential(&e_v, &te[0])) return false; //the te are modified inside approximate_potential()
      if(order >= 2) {
        e_v.x += step; 
        if(! approximate_potential(&e_v, &te2[0])) return false;
      }
      if(order >= 3) {
        e_v.x += step; 
        if(! approximate_potential(&e_v, &te3[0])) return false;
      }
      e_v = *R;
      e_v.x -= step; //find potential backward one step
      if(! approximate_potential(&e_v, &te[1])) return false;
      if(order >= 2) {
        e_v.x -= step;
        if(! approximate_potential(&e_v, &te2[1])) return false;
      }
      if(order >= 3) {
        e_v.x -= step;
        if(! approximate_potential(&e_v, &te3[1])) return false;
      }
      //cout << " _Scale in approx_force() = " << scale << endl; // New 3/22
      if(order == 1) dU[0] = scale * (te[1] - te[0]) / twostep;
      if(order == 2) dU[0] = scale * ((8. * (te[1] - te[0])) - (te2[1] - te2[0])) / twelvestep;
      if(order == 3) dU[0] = scale * ((45. * (te[1] - te[0])) - (9. * (te2[1] - te2[0])) + (te3[1] - te3[0])) / sixtystep;

      /* y axis */
      e_v = *R;
      e_v.y += step; 
      if(! approximate_potential(&e_v, &te[0])) return false;
      if(order >= 2) {
        e_v.y += step; 
        if(! approximate_potential(&e_v, &te2[0])) return false;
      }
      if(order >= 3) {
        e_v.y += step; 
        if(! approximate_potential(&e_v, &te3[0])) return false;
      }
      e_v = *R;
      e_v.y -= step;
      if(! approximate_potential(&e_v, &te[1])) return false;
      if(order >= 2) {
        e_v.y -= step;
        if(! approximate_potential(&e_v, &te2[1])) return false;
      }
      if(order >= 3) {
        e_v.y -= step;
        if(! approximate_potential(&e_v, &te3[1])) return false;
      }
      if(order == 1) dU[1] = scale * (te[1] - te[0]) / twostep;
      if(order == 2) dU[1] = scale * ((8. * (te[1] - te[0])) - (te2[1] - te2[0])) / twelvestep;
      if(order == 3) dU[1] = scale * ((45. * (te[1] - te[0])) - (9. * (te2[1] - te2[0])) + (te3[1] - te3[0])) / sixtystep;

      /* z axis */
      e_v = *R;
      e_v.z += step; 
      if(! approximate_potential(&e_v, &te[0])) return false;
      if(order >= 2) {
        e_v.z += step; 
        if(! approximate_potential(&e_v, &te2[0])) return false;
      }
      if(order >= 3) {
        e_v.z += step; 
        if(! approximate_potential(&e_v, &te3[0])) return false;
      }
      e_v = *R;
      e_v.z -= step;
      if(! approximate_potential(&e_v, &te[1])) return false;
      if(order >= 2) {
        e_v.z -= step;
        if(! approximate_potential(&e_v, &te2[1])) return false;
      }
      if(order >= 3) {
        e_v.z -= step;
        if(! approximate_potential(&e_v, &te3[1])) return false;
      }
      if(order == 1) dU[2] = scale * (te[1] - te[0]) / twostep;
      if(order == 2) dU[2] = scale * ((8. * (te[1] - te[0])) - (te2[1] - te2[0])) / twelvestep;
      if(order == 3) dU[2] = scale * ((45. * (te[1] - te[0])) - (9. * (te2[1] - te2[0])) + (te3[1] - te3[0])) / sixtystep;

      //cout << "U = " << te[1] << " kcal/mol\n";

      if(e) *e += E;

      F->x += dU[0];
      F->y += dU[1];
      F->z += dU[2];

      return true;
    }
   
    void translate(double dx, double dy, double dz) {
      origin[0] += dx;
      origin[1] += dy;
      origin[2] += dz;
    }

};




#endif
