/*--------------------------------------------------------------------------
  input:
    angle       - angle of rotation [rad]

  output:
    rotmat      - vector result
--------------------------------------------------------------------------*/
#include <cmath>
#include "R_x.h"
#include "Sat_Const.h"

void R_x(double angle, double rotmat[3][3]){

    double C = cos(angle);
    double S = sin(angle);

    rotmat[0][0] = 1.0;  rotmat[0][1] =    0.0;  rotmat[0][2] = 0.0;
    rotmat[1][0] = 0.0;  rotmat[1][1] =      C;  rotmat[1][2] =   S;
    rotmat[2][0] = 0.0;  rotmat[2][1] = -1.0*S;  rotmat[2][2] =   C;
}