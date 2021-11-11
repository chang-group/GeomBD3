#include "Main.h"
#include "Strings.h"
#include "../Gridder/Gridder.h"
#include "Timer.h"


void usage() {
  printf("Usage: ProbDX -r [Receptor.PDBQE] -t [Trajectory.PDB] -o [Output.DX] (Optional: -p GRID_PADDING(=10A) -s GRID_SPACING(=1.000A))\n");
}


int main(int argc, char **argv) {
  string datfn, trjfn, recfn, outfn, fldfn, stoken, token;
  double grid_resolution = 1.00;
  double padding = 10., padding_sqr;

  // receptor pdbqt
  if(!getInputWithFlag(argc, argv, 'r', &recfn)) { usage(); return -1; }
  // trj pdbqt
  if(!getInputWithFlag(argc, argv, 't', &trjfn)) { usage(); return -1; }
  // out pdbqt
  if(!getInputWithFlag(argc, argv, 'o', &outfn)) { usage(); return -1; }
  // number of processors/threads
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
  }
  // grid padding
  if(getInputWithFlag(argc, argv, 'p', &stoken)) {
    padding = stringToDouble(stoken);
  }
  padding_sqr = padding * padding;
  // grid spacing
  if(getInputWithFlag(argc, argv, 's', &stoken)) {
    grid_resolution = stringToDouble(stoken);
  }

  // Load receptor file
  cout << "> Loading receptor PDBQE..." << endl;
  Molecule_PQRE *rec = new Molecule_PQRE(recfn, NULL);
  cout << "> Receptor center: " << rec->center.x << ", " << rec->center.y << ", " << rec->center.z << endl;
  cout << "> Receptor minimum: " << rec->min.x << ", " << rec->min.y << ", " << rec->min.z << endl;
  cout << "> Receptor maximum coordinates: " << rec->max.x << ", " << rec->max.y << ", " << rec->max.z << endl;

  // Calculate grid geometries
  vertex origin = { rec->min.x - padding, rec->min.y - padding, rec->min.z - padding };
  vertex dimensions = { (rec->max.x + padding) - origin.x, (rec->max.y + padding) - origin.y, (rec->max.z + padding) - origin.z };
  int Npoints[3] = { dimensions.x / grid_resolution, dimensions.y / grid_resolution, dimensions.z / grid_resolution };
  int Ntotal = Npoints[0] * Npoints[1] * Npoints[2];

  //// Data for all grids
  long ***data = (long***)calloc(Npoints[0], sizeof(long**));
  for(int i=0; i < Npoints[0]; i++) {
    data[i] = (long**)calloc(Npoints[1], sizeof(long*));
    for(int j=0; j < Npoints[1]; j++) {
      data[i][j] = (long*)calloc(Npoints[2], sizeof(long));
    }
  }

  // Start a timer and start our calculations
  string line;
  ifstream trj;
  trj.open(trjfn.c_str(), ifstream::in);
  while(getline(trj, line)) {
    if(starts_with(&line, "ATOM")) {
      double x = stringToDouble(line.substr(30, 8));
      double y = stringToDouble(line.substr(38, 8));
      double z = stringToDouble(line.substr(46, 8));
      int Grec[3];
      if(coordinateToGrid(x, y, z, &Grec[0], &Grec[1], &Grec[2], &origin, grid_resolution)) {
        if(Grec[0] < 0 or Grec[0] >= Npoints[0] or
           Grec[1] < 0 or Grec[1] >= Npoints[1] or
           Grec[2] < 0 or Grec[2] >= Npoints[2]) {
          //...
        } else {
          data[Grec[0]][Grec[1]][Grec[2]] += 1;
        }
      }
    }
  }
  trj.close();
  cout << "> Done computing grid." << endl;

  FILE *fdo;
  fdo = fopen(outfn.c_str(), "w");
  fprintf(fdo, "object 1 class gridpositions counts %d %d %d\n", Npoints[0], Npoints[1], Npoints[2]);
  fprintf(fdo, "origin %12.6e %12.6e %12.6e\n", origin.x, origin.y, origin.z);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", grid_resolution, 0., 0.);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", 0., grid_resolution, 0.);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", 0., 0., grid_resolution);
  fprintf(fdo, "object 2 class gridconnections counts %d %d %d\n", Npoints[0], Npoints[1], Npoints[2]);
  fprintf(fdo, "object 3 class array type double rank 0 items %d data follows\n", Ntotal);
  int i[3] = { 0, 0, 0 };
  for(int it=0; it < Ntotal; it++) {
    fprintf(fdo, "%12.6e ", (float)data[i[0]][i[1]][i[2]]);
    if((it+1) % 3 == 0) fprintf(fdo, "\n");
    i[2]++;
    if(i[2] >= Npoints[2]) { i[2] = 0; i[1]++; }
    if(i[1] >= Npoints[1]) { i[1] = 0; i[0]++; }
  }
  if(Ntotal % 3 != 0) fprintf(fdo, "\n");
  fprintf(fdo, "attribute \"dep\" string \"positions\"\n");
  fprintf(fdo, "object \"regular positions regular connections\" class field\n");
  fprintf(fdo, "component \"positions\" value 1\n");
  fprintf(fdo, "component \"connections\" value 2\n");
  fprintf(fdo, "component \"data\" value 3\n");
  fclose(fdo);

  
  cout << "> Done writing OpenDX map." << endl;
  return EXIT_SUCCESS;
}
