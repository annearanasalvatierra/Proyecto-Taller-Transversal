/*--------------------------------------------------------------------------

 PrecMatrix.m

 Purpose:
   Precession transformation of equatorial coordinates

 Input:
   Mjd_1     Epoch given (Modified Julian Date TT)
   MjD_2     Epoch to precess to (Modified Julian Date TT)
 
 Output:
   PrecMat   Precession transformation matrix

 Last modified:   2015/08/12   M. Mahooti
 
--------------------------------------------------------------------------*/
#include <cmath>
#include "PrecMatrix.h"
#include "Sat_Const.h"
#include "Matriz.h"
#include "R_z.h"
#include "R_y.h"

void PrecMatrix (double Mjd_1, double Mjd_2, double PrecMat[3][3]){

    double T  = (Mjd_1-MJD_J2000)/36525.0;
    double dT = (Mjd_2-Mjd_1)/36525.0;

    // Precession angles
    double zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Arcs;
    double z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Arcs;
    double theta =  ( (2004.3109-(0.85330+0.000217*T)*T)-((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Arcs;

    // Precession matrix
    double matriz1[3][3] = {0};
    R_z(-z, matriz1);
    double matriz2[3][3] = {0};
    R_y(theta, matriz2);
    double matriz3[3][3] = {0};
    R_z(-zeta, matriz3);

    TripleProducto(matriz1, matriz2, matriz3, PrecMat);
}