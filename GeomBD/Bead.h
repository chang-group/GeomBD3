#ifndef Bead_h_
#define Bead_h_

#include "Main.h"


 
class Bead {
  private:
    vertex _R;
    vertex _F;

  public:
    vertex R;
    vertex F;
    double q;
    double r;
    double m;
    double S;
    string type;

  public:
    Bead() {
      R.x = 0.;
      R.y = 0.;
      R.z = 0.;
      F.x = 0.;
      F.y = 0.;
      F.z = 0.;
      q = 0.;
      r = 0.;
      m = 0.;
      S = 0.;
    }


    ~Bead() {
    }

    //Save data in public members to private members
    void save() {
      _R.x = R.x;
      _R.y = R.y;
      _R.z = R.z;
      _F.x = F.x;
      _F.y = F.y;
      _F.z = F.z;
    }

    //Restore public members data using data stored in private members
    void restore() {
      R.x = _R.x;
      R.y = _R.y;
      R.z = _R.z;
      F.x = _F.x;
      F.y = _F.y;
      F.z = _F.z;
    }

    //Translate vertex R by cartesian coordinates x,y,z
    void translate(double dx, double dy, double dz) {
      R.x += dx;
      R.y += dy;
      R.z += dz;
    }

    //Rotate vertex R. dax = difference in angle with x axis
    void rotate(double dax, double day, double daz) {
      double rm[3][3];
      double cosa, sina;
      double cosb, sinb;
      double cosc, sinc;

      double v[3] = { R.x, R.y, R.z };
      double p[3] = { 0., 0., 0. };

      cosa = cos(dax);
      sina = sin(dax);
      cosb = cos(day);
      sinb = sin(day);
      cosc = cos(daz);
      sinc = sin(daz);

      rm[0][0] = cosb*cosc;
      rm[0][1] = cosc*sina*sinb - cosa*sinc;
      rm[0][2] = cosa*cosc*sinb + sina*sinc;

      rm[1][0] = cosb*sinc;
      rm[1][1] = cosa*cosc + sina*sinb*sinc;
      rm[1][2] = -cosc*sina + cosa*sinb*sinc;

      rm[2][0] = -sinb;
      rm[2][1] = cosb*sina;
      rm[2][2] = cosa*cosb;

      for(int n=0; n < 3; n++)
        for(int m=0; m < 3; m++)
          p[n] += rm[n][m] * v[m];

      for(int n=0; n < 3; n++) {
        v[n] = p[n];
      }

      R.x = v[0];
      R.y = v[1];
      R.z = v[2];
    }


    void rotateAbout(double ox, double oy, double oz, double dax, double day, double daz) {
      translate(-ox, -oy, -oz);
      rotate(dax, day, daz);
      translate(ox, oy, oz);
    }

};


#endif
