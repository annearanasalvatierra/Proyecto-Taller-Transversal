/*------------------------------------------------------------------------------

 VarEqn.m

 Purpose:
   Computes the variational equations, i.e. the derivative of the state vector
   and the state transition matrix

 Input:
   x           Time since epoch in [s]
   yPhi        (6+36)-dim vector comprising the state vector (y) and the
               state transition matrix (Phi) in column wise storage order

 Output:
   yPhip       Derivative of yPhi
 
 Last modified:   2015/08/12   M. Mahooti

------------------------------------------------------------------------------*/
#include "VarEqn.h"
#include "Sat_Const.h"
#include "IERS.h"
#include "timediff.h"
#include "Matriz.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"

extern double **eopdata;

typedef struct{
	double Mjd_TT;
	int n;
	int m;
	double Mjd_UTC;
} Param;

extern Param AuxParam;

void VarEqn(double x, double yPhi[42], double yPhip[42]){

    double UT1_UTC, TAI_UTC, x_pole, y_pole;
    IERS (eopdata, AuxParam.Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    // Transformation matrix
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

    double E[3][3];
    TripleProducto(PoleMat, GHAMat, T, E);

    // State vector components
    double r[3];
    r[0] = yPhi[0];
    r[1] = yPhi[1];
    r[2] = yPhi[2];

    double v[3];
    v[0] = yPhi[3];
    v[1] = yPhi[4];
    v[2] = yPhi[5];

    double Phi[6][6] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};

    // State transition matrix
    for (int j=0; j<6; j++){ 
        int startIndex = 6 * j + 1; 
        int endIndex = 6 * j + 6; 

        for (int i = 0; i < 6; i++) {
            Phi[i][j] = yPhi[startIndex + i];
        }
    }

    // Acceleration and gradient
    double a[3];
    AccelHarmonic ( r, E, AuxParam.n, AuxParam.m, a);
    double G[3][3];
    G_AccelHarmonic ( r, E, AuxParam.n, AuxParam.m, G);

    // Time derivative of state transition matrix
    for (int i=0; i<42; i++){
        yPhip[i] = 0;
    }
    double dfdy[6][6] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            dfdy[i][j] = 0.0;                 // dv/dr(i,j)
            dfdy[i+3][j] = G[i][j];            // da/dr(i,j)
            if ( i==j ){
                dfdy[i][j+3] = 1;
            }else{
                dfdy[i][j+3] = 0;             // dv/dv(i,j)
            }
            dfdy[i+3][j+3] = 0.0;             // da/dv(i,j)
        }
    }

    double Phip[6][6] = {0};
    for (int z = 0; z < 6; z++) {
         for (int y = 0; y < 6; y++) {
                Phip[z][y] += dfdy[z][y] * Phi[z][y];
         }
    }

    // Derivative of combined state vector and state transition matrix
    for (int ii=0; ii<3; ii++){
        yPhip[ii]   = v[ii];                 // dr/dt(ii)
        yPhip[ii+3] = a[ii];                 // dv/dt(ii)
    }

    for (int iii=0; iii<6; iii++){
        for (int jj=0; jj<6; jj++){
            yPhip[6*jj+iii] = Phip[iii][jj];       // dPhi/dt(iii,jj)
        }
    }
}