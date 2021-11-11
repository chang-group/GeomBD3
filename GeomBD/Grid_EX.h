#ifndef _Grid_EX_h_
#define _Grid_EX_h_

#include "Main.h"
#include "Strings.h"

 
class Grid_EX {
  friend class Model;
  protected:
    bool header = true,
         body = false,
         footer = false;

    int N[3];
    int Nt;
    double origin[3];
    double delta;
    bool ***data;

  public:
    string type;

  public:
    Grid_EX(string bpm_filename) {
      ifstream fd(bpm_filename.c_str(), ios::in | ios::binary);

      if(!fd.is_open()) {
        cout << "! Failed to load grid file. File not found." << endl;
        exit(-1);
      }

      fd.read((char*)&origin, sizeof(double) * 3); //Read origin
      fd.read((char*)&N, sizeof(int) * 3);         //Read # of points in each dimension
      fd.read((char*)&delta, sizeof(double));      //Read spacing
      Nt = N[0] * N[1] * N[2];                     //Calc total points
//Allocate memory for grid by looping over all x y z points
      data = (bool***)calloc(N[0], sizeof(bool**));
      for(int nx=0; nx < N[0]; nx++) {
        data[nx] = (bool**)calloc(N[1], sizeof(bool*));              
        if(!data[nx]) {
          cout << "failed allocation\n";
          exit(EXIT_FAILURE);
        }
        for(int ny=0; ny < N[1]; ny++) {
          data[nx][ny] = (bool*)calloc(N[2], sizeof(bool));              
          if(!data[nx][ny]) {
            cout << "failed allocation\n";
            exit(EXIT_FAILURE);
          }
        }
      }

//Read all grid point values
      for(int nx=0; nx < N[0]; nx++) {
        for(int ny=0; ny < N[1]; ny++) {
          for(int nz=0; nz < N[2]; nz++) {
            fd.read((char*)&data[nx][ny][nz], sizeof(bool));
          }
        }
      }

    }

    virtual ~Grid_EX() {
      for(int nx=0; nx < N[0]; nx++) {
        for(int ny=0; ny < N[1]; ny++) {
          free(data[nx][ny]);
        }
        free(data[nx]);
      }
      free(data);
    }

    bool coordinateToGrid(vertex *R, int *G) {
      G[0] = (int)floor(((R->x-origin[0])/delta) + 0.5); //G = how many deltas (grid spacings) ligand is away from origin
      if(G[0] < 0 or G[0] >= N[0]) { //If ligand past grid on either side, return false. If G[0] < 0, we're before origin,
                                     //if >= N[0] (number of grid points), we're past edge
        return false;
      }
      G[1] = (int)floor(((R->y-origin[1])/delta) + 0.5);
      if(G[1] < 0 or G[1] >= N[1]) {
        return false;
      }
      G[2] = (int)floor(((R->z-origin[2])/delta) + 0.5);
      if(G[2] < 0 or G[2] >= N[2]) {
        return false;
      }

      return true;
    }

    bool onGrid(vertex *R) {
      int grid[3];
      return coordinateToGrid(R, (int*)&grid);
    }

    short value(vertex *R) {
      int grid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return 0;
      return (data[grid[0]][grid[1]][grid[2]] ? 1 : 0); //If data array holds value at point where ligand is, then return 1, else return 0?
    }

    void translate(double dx, double dy, double dz) {
      origin[0] += dx;
      origin[1] += dy;
      origin[2] += dz;
    }

};

#endif
