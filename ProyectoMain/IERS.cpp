/*--------------------------------------------------------------------------

 IERS.m

 Purpose:
   Management of IERS time and polar motion data
  
 Last modified:   2015/08/12   M. Mahooti
 
--------------------------------------------------------------------------*/
#include <cmath>
#include "Sat_Const.h"
#include "IERS.h"

void IERS (double **eopdata, double Mjd_TT, double& UT1_UTC, double& TAI_UTC, double& x_pole, double& y_pole){

    double Arcs = 3600.0*180.0/pi;  // Arcseconds per radian

    double mj = round(Mjd_TT);
    int nop = 13;

    double eopColumna[13];

    for (int i=0; i < nop; i++){
        if (mj == eopdata[3][i]){
            for (int j=0; j<nop; j++){
                eopColumna[j]=eopdata[j][i];
            }
        }      
    }

    // Setting of IERS Earth rotation parameters
    // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    UT1_UTC = eopColumna[6];      // UT1-UTC time difference [s]
    TAI_UTC = eopColumna[12];     // TAI-UTC time difference [s]
    x_pole  = eopColumna[4]/Arcs; // Pole coordinate [rad]
    y_pole  = eopColumna[5]/Arcs; // Pole coordinate [rad]
}