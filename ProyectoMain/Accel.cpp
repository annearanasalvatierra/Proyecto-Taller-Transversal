/*--------------------------------------------------------------------------
 Accel.m

 Purpose:
   Computes the acceleration of an Earth orbiting satellite due to 
    - the Earth's harmonic gravity field, 
    - the gravitational perturbations of the Sun and Moon
    - the solar radiation pressure and
    - the atmospheric drag

 Inputs:
   Mjd_TT      Terrestrial Time (Modified Julian Date)
   Y           Satellite state vector in the ICRF/EME2000 system

 Output:
   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system

 Last modified:   2015/08/12   M. Mahooti
 
--------------------------------------------------------------------------*/
#include <cmath>
#include "Accel.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "Sat_Const.h"
#include "Matriz.h"

extern double **eopdata;

typedef struct{
	double Mjd_TT;
	int n;
	int m;
	double Mjd_UTC;
} Param;

extern Param AuxParam;

void Accel(double x, double Y[6], double dY[6]) {

    double UT1_UTC, TAI_UTC, x_pole, y_pole;
    IERS(eopdata, AuxParam.Mjd_TT + x/86400, UT1_UTC, TAI_UTC, x_pole, y_pole);

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    double Mjd_UT1 = AuxParam.Mjd_TT + x/86400 + (UT1_UTC-TT_UTC)/86400.0;

    double P[3][3] = {0}; 
    PrecMatrix(MJD_J2000, AuxParam.Mjd_TT + x/86400, P);

    double N[3][3] = {0};
    NutMatrix(AuxParam.Mjd_TT + x/86400, N);

    double T[3][3] = {0};
    Producto3X3y3X3(N, P, T);

    double PoleMat[3][3] = {0};
    PoleMatrix(x_pole, y_pole, PoleMat);

    double GHAMat[3][3] = {0};
    GHAMatrix(Mjd_UT1, GHAMat);

    double E[3][3] = {0};
    TripleProducto(PoleMat, GHAMat, T, E);

    double a[3] = {0};
    double YTres[3];
    YTres[0] = Y[0];
    YTres[1] = Y[1];
    YTres[2] = Y[2];
    AccelHarmonic(YTres, E, AuxParam.n, AuxParam.m, a);

    dY[0] = Y[3];
    dY[1] = Y[4];
    dY[2] = Y[5];
    dY[3] = E[0][0];
    dY[4] = E[1][1];
    dY[5] = E[2][2];
}