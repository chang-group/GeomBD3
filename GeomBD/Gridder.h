#ifndef _Gridder_h_
#define _Gridder_h_
#include "Main.h"
#include "Strings.h"
#include "Timer.h"
using namespace std;


inline bool getInputWithFlag(int argc, char **argv, char flag, string *value) {
  int  opt;                                                                                                                                                                                                     
  bool bopt = false;
  char gopt[2] = { flag, ':' };

  for(int i=1; i < argc; i++) {
    if(argv[i][0] == '-' && argv[i][1] == flag) {
      if(argv[i+1][0] != '-') {
        *value = argv[i+1];
        bopt = true;
        i++;
        break;
      }
    }
  }

  return bopt;
}


inline bool coordinateToGrid(double x, double y, double z, int *Gx, int *Gy, int *Gz, vertex *origin, double delta) {
  *Gx = (int)floor(((x-origin->x)/delta) + 0.5);
  *Gy = (int)floor(((y-origin->y)/delta) + 0.5);
  *Gz = (int)floor(((z-origin->z)/delta) + 0.5);
  return true;
}


inline bool gridToCoordinate(int Gx, int Gy, int Gz, double *x, double *y, double *z, vertex *origin, double delta) {
  *x = (((double)Gx) * delta) + origin->x;
  *y = (((double)Gy) * delta) + origin->y;
  *z = (((double)Gz) * delta) + origin->z;
  return true;
}


struct LJ_Pair_Parameter {
  double A;
  double B;
  double e;
  double s;
};
typedef vector<LJ_Pair_Parameter> Vec;
typedef vector<Vec> Pair_Map;


class GBD3Parameters {
  public:
    vector<string> types;
    vector<double> Rii;
    vector<double> epsii;
    vector<double> vol;
    vector<double> solpar;

    Pair_Map lj_map;
    Pair_Map slj_map;

  public:
    GBD3Parameters(string filename) {
      string line, token;
      ifstream fd;

      fd.open(filename.c_str(), ifstream::in);
      while(getline(fd, line)) {
        if(line[0] == '#') continue;
        parseNextValue(&line, &token);

        if(token.length() > 2) {
          if(token == "atom_par") {
            parseNextValue(&line, &token);
            types.push_back(token);
            parseNextValue(&line, &token);
            Rii.push_back(stringToDouble(token));
            parseNextValue(&line, &token);
            epsii.push_back(stringToDouble(token));
            parseNextValue(&line, &token);
            vol.push_back(stringToDouble(token));
            parseNextValue(&line, &token);
            solpar.push_back(stringToDouble(token));
          }
        }
      }
//NOTES - creates a Pair_Map with a number of Vec vectors equal to the number of elements in the types vector
      lj_map = Pair_Map(types.size());
      slj_map = Pair_Map(types.size());

      for(int i=0; i < types.size(); i++) {
        for(int j=0; j < i+1; j++) {
          LJ_Pair_Parameter parm, sparm;

          double rij = 0.5 * (Rii[i] + Rii[j]);
          double eij = sqrt(epsii[i] * epsii[j]);
          parm = {
            /*A*/4. * eij * pow(rij, 12),
            /*B*/4. * eij * pow(rij, 6),
            eij,
            rij
          };
          lj_map[i].push_back(parm);

          sparm = {
            /*A*/eij * pow(rij, 12),
            /*B*/2. * eij * pow(rij, 6),
            eij,
            rij
          };
          slj_map[i].push_back(sparm);
        }
      }
    }

    int index_for_type(string type) {
      for(int i=0; i < types.size(); i++) {
        if(types[i] == type) return i;
      }
      cout << "* Warning: GBD3Parameters::index_for_type returning -1" << endl;
      return -1;
    }

};


class Molecule_PQRE {
  public:
    vector<vertex> coordinates;
    vertex center;
    vertex max;
    vertex min;
    vector<int> types;
    set<string> types_set;
    vector<double> charges;
    vector<double> radii;
    GBD3Parameters *params;

    Molecule_PQRE(string filename, GBD3Parameters *param) {
      string line, token;
      ifstream fd;

      params = param;
      center.x = center.y = center.z = 0.;
      max.x = max.y = max.z = -1e9;
      min.x = min.y = min.z =  1e9;

      fd.open(filename.c_str(), ifstream::in);
      while(getline(fd, line)) {
        if(line[0] == '#') continue;
        if(line.substr(0, 4) == "ATOM" or line.substr(0, 6) == "HETATM") {
          // coordinates
          vertex R;
          R.x = stringToDouble(line.substr(30, 10));
          R.y = stringToDouble(line.substr(40, 10));
          R.z = stringToDouble(line.substr(50, 10));
          coordinates.push_back(R);
          center.x += R.x;
          center.y += R.y;
          center.z += R.z;
          if(R.x > max.x) max.x = R.x;
          if(R.y > max.y) max.y = R.y;
          if(R.z > max.z) max.z = R.z;
          if(R.x < min.x) min.x = R.x;
          if(R.y < min.y) min.y = R.y;
          if(R.z < min.z) min.z = R.z;
          // type
          string t = line.substr(12, 4);
          string tt = trim(t);
          if(params)
            types.push_back(params->index_for_type(tt));
          types_set.insert(tt);
          // charge
          double q = stringToDouble(line.substr(60, 8));
          charges.push_back(q);
          // radius
          double r = stringToDouble(line.substr(68, 8));
          radii.push_back(r);
        }
      }

      center.x /= charges.size();
      center.y /= charges.size();
      center.z /= charges.size();
    }

    void print_types() {
      cout << "> Atom types in molecule:";
      set<string>::iterator it;
      for(it=types_set.begin(); it!=types_set.end(); ++it) {
        cout << ' ' << *it;
      }
      cout << endl;
    }

    bool check_types() {
      if(params) {
        set<string>::iterator it;
        for(it=types_set.begin(); it!=types_set.end(); ++it) {
          bool found = false;
          for(int j=0; j < params->types.size(); j++) {
            if(params->types[j] == *it) found = true;
          }
          if(!found) {
            cout << "! Error: Ligand atom type " << *it << " not found in the parameters." << endl;
            return false;
          }
        }
      } else {
        cout << "! Can't check for molecule atom types for exitence in parameter set, as no parameter set was specified." << endl;
      }
      return true;
    }
};



class Map_Potential {
  public:
    vertex origin;
    vertex dimensions;
    int Npoints[3];
    long long int Ntotal;

    Molecule_PQRE *receptor;
    double padding;
    double spacing;

    double* data_t;
    int type_t;
    ofstream *bpm_t;
    string filename;

    Map_Potential(string bpmfn, int atom_type_index, Molecule_PQRE *rec, double Arg_GridSpacing, double Arg_Padding) {
      origin = { rec->min.x - Arg_Padding, rec->min.y - Arg_Padding, rec->min.z - Arg_Padding };
      dimensions = { (rec->max.x + Arg_Padding) - origin.x, (rec->max.y + Arg_Padding) - origin.y, (rec->max.z + Arg_Padding) - origin.z };
      Npoints[0] = dimensions.x / Arg_GridSpacing;
      Npoints[1] = dimensions.y / Arg_GridSpacing;
      Npoints[2] = dimensions.z / Arg_GridSpacing;
      Ntotal = Npoints[0] * Npoints[1] * Npoints[2];

      receptor = rec;
      padding = Arg_Padding;
      spacing = Arg_GridSpacing;

      data_t = (double*)calloc(Npoints[2], sizeof(double));
      if(!data_t) { cout << "! Error: Could not allocate memory for grid calculation." << endl; exit(EXIT_FAILURE); }
      filename = bpmfn;
      type_t = atom_type_index;

      // open file
      bpm_t = new ofstream(filename, ios::out | ios::binary);
      // write header
      bpm_t->write((char*)&origin, sizeof(vertex));
      bpm_t->write((char*)&Npoints, sizeof(int)*3);
      bpm_t->write((char*)&Arg_GridSpacing, sizeof(double));
    }

    ~Map_Potential() {
      bpm_t->close();
      delete bpm_t;
      free(data_t);
    }
};


class Map_Exclusion {
  public:
    bool ***data;
    vertex origin;
    vertex dimensions;
    int Npoints[3];
    long long int Ntotal;

    Molecule_PQRE *receptor;
    double padding;
    double spacing;
    double scaling;
    bool include_padding;

    Map_Exclusion(Molecule_PQRE *rec, double Arg_GridSpacing, double Arg_Padding, double Arg_Scaling=1.0, bool Use_Padding=true) {
      origin = { rec->min.x - Arg_Padding, rec->min.y - Arg_Padding, rec->min.z - Arg_Padding };
      dimensions = { (rec->max.x + Arg_Padding) - origin.x, (rec->max.y + Arg_Padding) - origin.y, (rec->max.z + Arg_Padding) - origin.z };
      Npoints[0] = dimensions.x / Arg_GridSpacing;
      Npoints[1] = dimensions.y / Arg_GridSpacing;
      Npoints[2] = dimensions.z / Arg_GridSpacing;
      Ntotal = Npoints[0] * Npoints[1] * Npoints[2];

      receptor = rec;
      padding = Arg_Padding;
      spacing = Arg_GridSpacing;
      scaling = Arg_Scaling;
      include_padding = Use_Padding;

      // allocate
      data = (bool***)calloc(Npoints[0], sizeof(bool**));
      for(int i=0; i < Npoints[0]; i++) {
        data[i] = (bool**)calloc(Npoints[1], sizeof(bool*));
        for(int j=0; j < Npoints[1]; j++) {
          data[i][j] = (bool*)calloc(Npoints[2], sizeof(bool));
        }
      }
    }


    Map_Exclusion(string bpmFN) {
      ifstream fd(bpmFN, ios::in | ios::binary);

      fd.read((char*)&origin, sizeof(double) * 3);
      fd.read((char*)&Npoints, sizeof(int) * 3);
      fd.read((char*)&spacing, sizeof(double));

      Ntotal = Npoints[0] * Npoints[1] * Npoints[2];

      receptor = NULL;
      padding = 0.;

      // allocate
      data = (bool***)calloc(Npoints[0], sizeof(bool**));
      for(int i=0; i < Npoints[0]; i++) {
        data[i] = (bool**)calloc(Npoints[1], sizeof(bool*));
        for(int j=0; j < Npoints[1]; j++) {
          data[i][j] = (bool*)calloc(Npoints[2], sizeof(bool));
        }
      }

      // read file
      for(int nx=0; nx < Npoints[0]; nx++) {
        for(int ny=0; ny < Npoints[1]; ny++) {
          for(int nz=0; nz < Npoints[2]; nz++) {
            fd.read((char*)&data[nx][ny][nz], sizeof(bool));
          }
        }
      }
    }


    void calculate() {
      if(receptor == NULL) {
        cout << "! 'calculate()' called on Map_Exclusion object with no defined receptor." << endl;
        return;
      }

      cout << "> Creating exclusion grid...";
      cout.flush();
      cilk_for(int i=0; i < receptor->coordinates.size(); i++) {
        vertex Rrec = receptor->coordinates[i];
        double radius = receptor->radii[i];
        int Grec[3];
        if(coordinateToGrid(Rrec.x, Rrec.y, Rrec.z, &Grec[0], &Grec[1], &Grec[2], &origin, spacing)) {
          double searchrd = 0.;
          int searchr;
          if(include_padding) {
            searchrd = (radius*scaling) + padding;
          } else {
            searchrd = radius*scaling;
          }
          searchr = (searchrd / spacing) + 1;
          for(int gx=(Grec[0]-searchr); gx <= Grec[0]+searchr; gx++) {
            if(gx < 0 or gx >= Npoints[0]) continue;
            for(int gy=(Grec[1]-searchr); gy <= Grec[1]+searchr; gy++) {
              if(gy < 0 or gy >= Npoints[1]) continue;
              for(int gz=(Grec[2]-searchr); gz <= Grec[2]+searchr; gz++) {
                if(gz < 0 or gz >= Npoints[2]) continue;

                double dx = (origin.x + (gx * spacing)) - Rrec.x;
                double dy = (origin.y + (gy * spacing)) - Rrec.y;
                double dz = (origin.z + (gz * spacing)) - Rrec.z;
                if(sqrt(dx*dx + dy*dy + dz*dz) <= searchrd) {
                  data[gx][gy][gz] = true;
                }
              }
            }
          }
        }
      }
      cout << " Done." << endl;
    }


    ~Map_Exclusion() {
      // deallocate
      for(int i=0; i < Npoints[0]; i++) {
        for(int j=0; j < Npoints[1]; j++) {
          free(data[i][j]);
        }
        free(data[i]);
      }
      free(data);
    }

    void write(string filename) {
      ofstream fo(filename);
      cout << "> Writing bxm..." << endl;
      fo.write((char*)&origin, sizeof(vertex));
      fo.write((char*)&Npoints, sizeof(int)*3);
      fo.write((char*)&spacing, sizeof(double));

      for(int i=0; i < Npoints[0]; i++) {
        for(int j=0; j < Npoints[1]; j++) {
          fo.write((char*)data[i][j], sizeof(bool) * Npoints[2]);
        }
      }
    }


    void write_dx(string filename) {
      printf("> Starting to write OpenDX potential map...\n");
      // output opendx data
      FILE *fdo;
      fdo = fopen(filename.c_str(), "w");
      fprintf(fdo, "object 1 class gridpositions counts %d %d %d\n", Npoints[0], Npoints[1], Npoints[2]);
      fprintf(fdo, "origin %12.6e %12.6e %12.6e\n", origin.x, origin.y, origin.z);
      fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", spacing, 0., 0.);
      fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", 0., spacing, 0.);
      fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", 0., 0., spacing);
      fprintf(fdo, "object 2 class gridconnections counts %d %d %d\n", Npoints[0], Npoints[1], Npoints[2]);
      fprintf(fdo, "object 3 class array type double rank 0 items %d data follows\n", Ntotal);
      int i[3] = { 0, 0, 0 };
      for(long long int it=0; it < Ntotal; it++) {
        fprintf(fdo, "%12.6e ", (double)data[i[0]][i[1]][i[2]]);
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
    }


    bool coordToGrid(double x, double y, double z, int *gx, int *gy, int *gz) {
      coordinateToGrid(x, y, z, gx, gy, gz, &origin, spacing);
      if(*gx < 0 or *gx > Npoints[0] or *gy < 0 or *gy > Npoints[1] or *gz < 0 or *gz > Npoints[2]) {
        return false;
      }
      return true;
    }


    bool gridToCoord(int gx, int gy, int gz, double *x, double *y, double *z) {
      if(gx < 0 or gx > Npoints[0] or gy < 0 or gy > Npoints[1] or gz < 0 or gz > Npoints[2]) {
        return false;
      }
      gridToCoordinate(gx, gy, gz, x, y, z, &origin, spacing);
      return true;
    }

    bool get_value(int gx, int gy, int gz) {
      if(gx < 0 or gy < 0 or gz < 0 or gx >= Npoints[0] or gy >= Npoints[1] or gz >= Npoints[2]) return false;
      return data[gx][gy][gz];
    }

};


#endif
