/*--------------------------------------------------------------------------

 AccelHarmonic.m

 Purpose:
   Computes the acceleration due to the harmonic gravity field of the 
   central body

 Inputs:
 Mjd_TT        Modified Julian Date of TT
   r           Satellite position vector in the inertial system
   E           Transformation matrix to body-fixed system
   n_max       Maximum degree
   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)

 Output:
   a           Acceleration (a=d^2r/dt^2)

 Last modified:   2015/08/12   M. Mahooti
 
--------------------------------------------------------------------------*/
#include "AccelHarmonic.h"
#include "Sat_Const.h"
#include "vector.h"
#include "Matriz.h"
#include "Legendre.h"
#include <cmath>

extern double **Cnm, **Snm;

void AccelHarmonic(double r[3], double E[3][3], int n_max, int m_max, double a[3]){

    double gm    = 398600.4415e9; // [m^3/s^2]; JGM3/EGM96
    double r_ref = 6378.1363e3;   // Radius Earth [m]; JGM3/EGM96

    // Body-fixed position 
    double r_bf[3];

    // Realizar el producto de la matriz E y del vector r
    for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
             r_bf[i] += E[i][k] * r[k];
         }
    }

    // Auxiliary quantities
    double d = norm(r_bf, 3);   // distance
    double latgc = asin(r_bf[2]/d);
    double lon = atan2(r_bf[1],r_bf[0]);

    double** pnm; 
    double** dpnm;
    Legendre(n_max, m_max, latgc, &pnm, &dpnm);

    double dUdr = 0;
    double dUdlatgc = 0;
    double dUdlon = 0;
 
    double q1 = 0;
    double q2 = 0;
    double q3 = 0; 

    double b1, b2, b3;

    for (int n=0; n<n_max; n++){
        b1 = (-gm/pow(d,2))*pow((r_ref/d),n*(n+1));
        b2 = (gm/d)*pow((r_ref/d),n);
        b3 = (gm/d)*pow((r_ref/d),n);
        for (int m=0; m<n; m++){
            q1 = q1 + pnm[n][m]*(Cnm[n][m]*cos(m*lon)+Snm[n][m]*sin(m*lon));
            q2 = q2 + dpnm[n][m]*(Cnm[n][m]*cos(m*lon)+Snm[n][m]*sin(m*lon));
            q3 = q3 + m*pnm[n][m]*(Snm[n][m]*cos(m*lon)-Cnm[n][m]*sin(m*lon));
        }
        dUdr     = dUdr     + q1*b1;
        dUdlatgc = dUdlatgc + q2*b2;
        dUdlon   = dUdlon   + q3*b3;

        q3 = 0; 
        q2 = q3; 
        q1 = q2;
    }

    // Body-fixed acceleration
    double r2xy = pow(r_bf[0],2)+pow(r_bf[1],2);

    double ax = (1/d*dUdr-r_bf[2]/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf[0]-(1/r2xy*dUdlon)*r_bf[1];
    double ay = (1/d*dUdr-r_bf[2]/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf[1]+(1/r2xy*dUdlon)*r_bf[0];
    double az =  1/d*dUdr*r_bf[2]+sqrt(r2xy)/pow(d,2)*dUdlatgc;

    double a_bf[3];
    a_bf[0] = ax;
    a_bf[1] = ay;
    a_bf[2] = az;

    // Inertial acceleration 
    double transpuesta[3][3];
    Transpuesta(E, transpuesta);

    // Realizar el producto de la matriz transpuesta y del vector a_bf
    for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
             a[i] += transpuesta[i][k] * a_bf[k];
         }
    }
}