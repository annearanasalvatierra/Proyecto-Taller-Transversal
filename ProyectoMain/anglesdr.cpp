/*---------------------------------------------------------------------------
 
  anglesdr.m

  this function solves the problem of orbit determination using three
  optical sightings.
 
  inputs:
    az1      - azimuth at t1               rad
    az2      - azimuth at t2               rad
    az3      - azimuth at t3               rad
    el1      - elevation at t1             rad
    el2      - elevation at t2             rad
    el3      - elevation at t3             rad
    Mjd1     - Modified julian date of t1
    Mjd2     - Modified julian date of t2
    Mjd3     - Modified julian date of t3
    rsite1   - ijk site1 position vector   m
    rsite2   - ijk site2 position vector   m
    rsite3   - ijk site3 position vector   m

  outputs:
    r        - ijk position vector at t2   m
    v        - ijk velocity vector at t2   m/s

 Last modified:   2015/08/12   M. Mahooti
 
---------------------------------------------------------------------------*/
#include "anglesdr.h"
#include "Sat_Const.h"
#include "Geodetic.h"
#include "LTC.h"
#include "Matriz.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "vector.h"
#include "doubler.h"
#include <cmath>

extern double **eopdata;

void anglesdr (double az1, double az2,double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, double rsite1[3], double rsite2[3], double rsite3[3], double r2[3], double v2[3]){

    double magr1in = 1.1*R_Earth;
    double magr2in = 1.11*R_Earth;
    char direct  = 'y';

    double tol    = 1e-8*R_Earth;
    double pctchg = 0.005;

    double t1 = (Mjd1 - Mjd2)*86400.0;
    double t3 = (Mjd3 - Mjd2)*86400.0;

    double los1[3];
    los1[0] = cos(el1)*sin(az1); 
    los1[1] = cos(el1)*cos(az1); 
    los1[2] = sin(el1);

    double los2[3];
    los2[0] = cos(el2)*sin(az2);
    los2[1] = cos(el2)*cos(az2); 
    los2[2] = sin(el2);

    double los3[3];
    los3[0] = cos(el3)*sin(az3);
    los3[1] = cos(el3)*cos(az3); 
    los3[2] = sin(el3);

    double lon1, lat1, h1;
    Geodetic(rsite1, lon1, lat1, h1);

    double lon2, lat2, h2;
    Geodetic(rsite2, lon2, lat2, h2);

    double lon3, lat3, h3;
    Geodetic(rsite3, lon3, lat3, h3);

    double M1[3][3] = {0};
    double M2[3][3] = {0};
    double M3[3][3] = {0};
    LTC(lon1, lat1, M1);
    LTC(lon2, lat2, M2);
    LTC(lon3, lat3, M3);

    // body-fixed system
    double transpuestaM1[3][3];
    Transpuesta(M1, transpuestaM1);

    // Realizar el producto de M1'*los1
    for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
              los1[i] += transpuestaM1[i][k] * los1[k];
         }  
    }

    // Realizar el producto de M1'*los2
    for (int ii = 0; ii < 3; ii++) {
         for (int kk = 0; kk < 3; kk++) {
              los2[ii] += transpuestaM1[ii][kk] * los2[kk];
         }  
    }

    // Realizar el producto de M1'*los3
    for (int iii = 0; iii < 3; iii++) {
         for (int kkk = 0; kkk < 3; kkk++) {
              los3[iii] += transpuestaM1[iii][kkk] * los3[kkk];
         }  
    }

    // mean of date system (J2000)
    double Mjd_UTC = Mjd1;
    double UT1_UTC, TAI_UTC, x_pole, y_pole;
    IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_TT = Mjd_UTC + TT_UTC/86400;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    double P[3][3] = {0};
    PrecMatrix(MJD_J2000, Mjd_TT, P);

    double N[3][3] = {0};
    NutMatrix(Mjd_TT, N);

    double T[3][3] = {0};
    Producto3X3y3X3(N, P, T);

    double PoleMat[3][3] = {0};
    PoleMatrix(x_pole, y_pole, PoleMat);

    double GHAMat[3][3] = {0};
    GHAMatrix(Mjd_UT1, GHAMat);

    double E[3][3] = {0};
    TripleProducto(PoleMat, GHAMat, T, E);

    double transpuestaE[3][3] = {0};
    Transpuesta(E, transpuestaE);

    // Realizar el producto de E'*los1
    for (int j = 0; j < 3; j++) {
         for (int l = 0; l < 3; l++) {
              los1[j] += transpuestaE[j][l] * los1[l];
         }  
    }

    // Realizar el producto de E'*rsite1
    for (int jj = 0; jj < 3; jj++) {
         for (int ll = 0; ll < 3; ll++) {
              rsite1[jj] += transpuestaE[jj][ll] * rsite1[ll];
         }  
    }

    Mjd_UTC = Mjd2;
    double UT1_UTC1, TAI_UTC1, x_pole1, y_pole1;
    IERS(eopdata, Mjd_UTC, UT1_UTC1, TAI_UTC1, x_pole1, y_pole1);

    double UT1_TAI1, UTC_GPS1, UT1_GPS1, TT_UTC1, GPS_UTC1;
    timediff(UT1_UTC1, TAI_UTC1, UT1_TAI1, UTC_GPS1, UT1_GPS1, TT_UTC1, GPS_UTC1);
    Mjd_TT = Mjd_UTC + TT_UTC1/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC1-TT_UTC1)/86400.0;

    double P1[3][3] = {0};
    PrecMatrix(MJD_J2000, Mjd_TT, P1);

    double N1[3][3] = {0};
    NutMatrix(Mjd_TT, N1);

    double T1[3][3] = {0};
    Producto3X3y3X3(N1, P1, T1);

    double PoleMat1[3][3] = {0};
    PoleMatrix(x_pole1, y_pole1, PoleMat1);

    double GHAMat1[3][3] = {0};
    GHAMatrix(Mjd_UT1, GHAMat1);

    double E1[3][3] = {0};
    TripleProducto(PoleMat1, GHAMat1, T1, E1);

    double transpuestaE1[3][3] = {0};
    Transpuesta(E1, transpuestaE1);

    // Realizar el producto de E1'*los2
    for (int aaa = 0; aaa < 3; aaa++) {
         for (int b = 0; b < 3; b++) {
              los2[aaa] += transpuestaE1[aaa][b] * los2[b];
         }  
    }

    // Realizar el producto de E1'*rsite2
    for (int aa = 0; aa < 3; aa++) {
         for (int bb = 0; bb < 3; bb++) {
              rsite2[aa] += transpuestaE1[aa][bb] * rsite2[bb];
         }  
    }

    Mjd_UTC = Mjd3;
    double UT1_UTC2, TAI_UTC2, x_pole2, y_pole2;
    IERS(eopdata, Mjd_UTC, UT1_UTC2, TAI_UTC2, x_pole2, y_pole2);
    double UT1_TAI2, UTC_GPS2, UT1_GPS2, TT_UTC2, GPS_UTC2;
    timediff(UT1_UTC2, TAI_UTC2, UT1_TAI2, UTC_GPS2, UT1_GPS2, TT_UTC2, GPS_UTC2);
    Mjd_TT = Mjd_UTC + TT_UTC2/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC2-TT_UTC2)/86400.0;

    double P2[3][3] = {0};
    PrecMatrix(MJD_J2000, Mjd_TT, P2);

    double N2[3][3] = {0};
    NutMatrix(Mjd_TT, N2);

    double T2[3][3] = {0};
    Producto3X3y3X3(N2, P2, T2);

    double PoleMat2[3][3] = {0};
    PoleMatrix(x_pole2, y_pole2, PoleMat2);

    double GHAMat2[3][3] = {0};
    GHAMatrix(Mjd_UT1, GHAMat2);

    double E2[3][3] = {0};
    TripleProducto(PoleMat2, GHAMat2, T2, E2);

    double transpuestaE2[3][3] = {0};
    Transpuesta(E2, transpuestaE2);

    // Realizar el producto de E2'*los3
    for (int c = 0; c < 3; c++) {
         for (int d = 0; d < 3; d++) {
              los3[c] += transpuestaE2[c][d] * los3[d];
         }  
    }

    // Realizar el producto de E2'*rsite3
    for (int cc = 0; cc < 3; cc++) {
         for (int dd = 0; dd < 3; dd++) {
              rsite3[cc] += transpuestaE2[cc][dd] * rsite3[dd];
         }  
    }

    double magr1old  = 99999999.9;
    double magr2old  = 99999999.9;
    double magrsite1 = norm(rsite1, 3);
    double magrsite2 = norm(rsite2, 3);
    double magrsite3 = norm(rsite3, 3);

    double cc1 = 2.0*dot(los1, rsite1, 3, 3);
    double cc2 = 2.0*dot(los2, rsite2, 3, 3);
    int ktr = 0;

    double f, g, magr1o, deltar1, pf1pr1, pf2pr1, magr2o, deltar2, pf1pr2, pf2pr2, delta, delta1, delta2;
    double f1, f2, q1, q2, q3, magr1, magr2, a, deltae32, f1delr1, f2delr1,f1delr2,f2delr2;
    double r3[3];

    while (fabs(magr1in-magr1old) > tol || fabs(magr2in-magr2old) > tol){
        ktr = ktr + 1; 
        r2[0] = 0;
        r2[1] = 0;
        r2[2] = 0;
        r3[0] = 0;
        r3[1] = 0;
        r3[2] = 0;
        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2,r3,f1,f2,q1,magr1,magr2,a,deltae32);
        f  = 1.0 - a/magr2*(1.0-cos(deltae32));
        g  = t3 - sqrt(pow(fabs(a),3)/GM_Earth)*(deltae32-sin(deltae32));
        v2[0] = (r3[0] - f*r2[0])/g;
        v2[1] = (r3[1] - f*r2[1])/g;
        v2[2] = (r3[2] - f*r2[2])/g;
    
        magr1o = magr1in;
        magr1in = (1.0+pctchg)*magr1in;
        deltar1 = pctchg*magr1in;
        r2[0] = 0;
        r2[1] = 0;
        r2[2] = 0;
        r3[0] = 0;
        r3[1] = 0;
        r3[2] = 0;
        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32);
        pf1pr1 = (f1delr1-f1)/deltar1;
        pf2pr1 = (f2delr1-f2)/deltar1;
 
        magr1in = magr1o;
        deltar1 = pctchg*magr1in;
        magr2o = magr2in;
        magr2in = (1.0+pctchg)*magr2in;
        deltar2 = pctchg*magr2in;

        r2[0] = 0;
        r2[1] = 0;
        r2[2] = 0;
        r3[0] = 0;
        r3[1] = 0;
        r3[2] = 0;
        doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2,r3,f1delr2,f2delr2,q3,magr1,magr2,a,deltae32);
        pf1pr2 = (f1delr2-f1)/deltar2;
        pf2pr2 = (f2delr2-f2)/deltar2;

        magr2in = magr2o;
        deltar2 = pctchg*magr2in;
    
        delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;

        delta1 = pf2pr2*f1 - pf1pr2*f2;
        delta2 = pf1pr1*f2 - pf2pr1*f1;
        //HAsta aqui bien

        deltar1 = -delta1/delta;
        deltar2 = -delta2/delta;
    
        magr1old = magr1in;
        magr2old = magr2in;
    
        magr1in = magr1in + deltar1;
        magr2in = magr2in + deltar2;
    }
    f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    g  = t3 - sqrt(pow(fabs(a),3)/GM_Earth)*(deltae32-sin(deltae32));
    v2[0] = (r3[0] - f*r2[0])/g;
    v2[1] = (r3[1] - f*r2[1])/g;
    v2[2] = (r3[2] - f*r2[2])/g;
}