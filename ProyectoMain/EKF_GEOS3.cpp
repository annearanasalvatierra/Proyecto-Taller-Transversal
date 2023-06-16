/*--------------------------------------------------------------------------

 Initial Orbit Determination using Double-R-Iteration and Extended Kalman
 Filter methods

 Last modified:   2015/08/12   M. Mahooti

 Refrences:
 
   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
   Applications", Springer Verlag, Heidelberg, 2000
   
   D. Vallado, "Fundamentals of Astrodynamics and Applications", 
   3rd Edition, 2007

   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003

--------------------------------------------------------------------------*/
#include "Sat_Const.h"
#include "Mjday.h"
#include "Position.h"
#include "anglesdr.h"
#include "IERS.h"
#include "timediff.h"
#include "LTC.h"
#include "gmst.h"
#include "R_z.h"
#include "MeasUpdate.h"
#include "vector.h"
#include "TimeUpdate.h"
#include "DEIntegAccel.h"
#include "DEIntegVarEqn.h"
#include "AzElPa.h"
#include "Globales.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

using namespace std;

double **Cnm, **Snm, **eopdata;
Param AuxParam;

int main()
{
    FILE *fp;
    int f, c;
    double aux1, aux2;

	int n_eqn;

    fp = fopen("./data/egm.txt", "r");
    if(fp == NULL){
		printf("Fail open GGM03S.txt file\n");
		exit(EXIT_FAILURE);
	}
    
    Cnm = (double **) malloc(362 * sizeof(double *));
    if(Cnm == NULL){
        printf("Cnm: memory not allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i <= 362; i++){
        Cnm[i] = (double *) malloc(362 * sizeof(double));
        if(Cnm[i] == NULL){
            printf("Cnm[i]: memory not allocated\n");
            exit(EXIT_FAILURE);
        }
    }
    
    Snm = (double **) malloc(362 * sizeof(double *));
    if(Snm == NULL){
        printf("Snm: memory not allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i <= 362; i++){
        Snm[i] = (double *) malloc(362 * sizeof(double));
        if(Snm[i] == NULL){
            printf("Snm[i]: memory not allocated\n");
            exit(EXIT_FAILURE);
        }
    }
    for(int n = 0; n <= 361; n++){
		for(int m = 0; m <= n; m++){
			fscanf(fp, "%d%d%lf%lf%lf%lf", &f, &c, &Cnm[n][m], &Snm[n][m], &aux1, &aux2);
		}
	}
	fclose(fp);
    //cout << Cnm[45][23] << endl;
    //cout << Snm[1][1] << endl;

    //read Earth orientation parameters
    //  ----------------------------------------------------------------------------------------------------
    // |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    // |(0h UTC)           "         "          s          s          "        "          "         "     s
    //  ----------------------------------------------------------------------------------------------------

    fp = fopen("./data/eop19620101.txt", "r");
    if(fp == NULL){
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    eopdata = (double **) malloc(13 * sizeof(double *));
    if(eopdata == NULL){
        printf("eopdata: memory not allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i <= 13; i++){
        eopdata[i] = (double *) malloc(19713 * sizeof(double));
        if(eopdata[i] == NULL){
            printf("eopdata[i]: memory not allocated\n");
            exit(EXIT_FAILURE);
        }
    }

    for(int i = 0; i < 19713; i++){
        fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &eopdata[0][i], &eopdata[1][i], &eopdata[2][i], &eopdata[3][i], &eopdata[4][i], &eopdata[5][i], &eopdata[6][i], &eopdata[7][i], &eopdata[8][i], &eopdata[9][i], &eopdata[10][i], &eopdata[11][i], &eopdata[12][i]);
    }
    fclose(fp);
    
    //cout << eopdata[0][0] << endl;

    int nobs = 18;
    double obs[18][4];

    // read observations
    fp = fopen("./data/GEOS3.txt", "r");
    if(fp == NULL){
        printf("Fail open GEOS3.txt file\n");
        exit(EXIT_FAILURE);
    }

    char line[55], y[5], mo[3], d[3], h[3], mi[3], sss[7], a[9], e[9], di[10];
    int YYY, M, D, hh, mm, az, el, Distt;
    double ss;
    
    for(int i = 0; i < 18; i++){
        fgets(line, sizeof(line)+2, fp);

        strncpy(y, &line[0], 4);
        y[4]='\0';
        YYY=atoi(y);

        strncpy(mo, &line[5], 2);
        mo[2]='\0';
        M=atoi(mo);

        strncpy(d, &line[8],2);
        d[2]='\0';
        D=atoi(d);

        strncpy(h, &line[12], 2);
        h[2]='\0';
        hh=atoi(h);

        strncpy(mi, &line[15], 2);
        mi[2]='\0';
        mm=atoi(mi);

        strncpy(sss, &line[18], 6);
        sss[6]='\0';
        ss=atof(sss);

        strncpy(a, &line[25], 8);
        a[8]='\0';
        az=atof(a);

        strncpy(e, &line[35], 8);
        e[8]='\0';
        el=atof(e);

        strncpy(di, &line[44], 9);
        di[9]='\0';
        Distt=atof(di);

        cout << YYY << "  " << M <<  "  " << D <<  "  " << hh <<  "  " << mm <<  "  " << ss <<  "  " << az <<  "  " << el <<  "  " << Distt << endl;

        obs[i][0] = Mjday(YYY,M,D,hh,mm,ss);
        obs[i][1] = Rad*az;
        obs[i][2] = Rad*el;
        obs[i][3] = 1e3*Distt;
    }
    fclose(fp);
    
    for(int i = 0; i <= 14; i++){
        free(eopdata[i]);
    }
    free(eopdata);


    double sigma_range = 92.5;     // [m]
    double sigma_az = 0.0224*Rad;  // [rad]
    double sigma_el = 0.0139*Rad;  // [rad]

    // Kaena Point station
    double lat = Rad*21.5748;     // [rad]
    double lon = Rad*(-158.2706); // [rad]
    double alt = 300.20;          // [m]

    double Rs[3];
    Position(lon, lat, alt, Rs);

    double Mjd1 = obs[0][0];
    double Mjd2 = obs[8][0];
    double Mjd3 = obs[17][0];

    double r2[3];
    double v2[3];
    anglesdr ( obs[0][1],obs[8][1],obs[17][1],obs[0][2],obs[8][2],obs[17][2],Mjd1,Mjd2,Mjd3,Rs,Rs,Rs,r2,v2 );

    double Y0_apr[6];
    Y0_apr[0] = r2[0];
    Y0_apr[1] = r2[1];
    Y0_apr[2] = r2[2];
    Y0_apr[3] = v2[0];
    Y0_apr[4] = v2[1];
    Y0_apr[5] = v2[2];

    double Mjd0 = Mjday(1995,1,29,2,38,00.0);

    double Mjd_UTC = obs[8][0];
    double UT1_UTC, TAI_UTC, x_pole, y_pole;
    IERS(eopdata,Mjd_UTC,UT1_UTC,TAI_UTC,x_pole,y_pole);

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    AuxParam.Mjd_TT  = Mjd_UTC + TT_UTC/86400; 
    AuxParam.n       = 10;
    AuxParam.m       = 10;

    n_eqn  = 6;

    DEIntegAccel(0,-(obs[8][0]-Mjd0)*86400.0, 1e-13, 1e-6, 6, Y0_apr);
    double Y[6][1];
    Y[0][0] = Y0_apr[0];
    Y[1][0] = Y0_apr[1];
    Y[2][0] = Y0_apr[2];
    Y[3][0] = Y0_apr[3];
    Y[4][0] = Y0_apr[4];
    Y[5][0] = Y0_apr[5];

    double P[6][6];
    for(int i=0; i<6; i++){
        for(int j=0; j<6; j++){
            P[i][j] = 0;
        }
    }
  
    for (int i=0; i<3; i++){
        P[i][i] = 1e8;
    }

    for (int i=3; i<6; i++){
        P[i][i] = 1e3;
    }

    double LT[3][3];
    LTC(lon,lat,LT);

    double yPhi[42];
    for(int j=0; j<42; j++){
        yPhi[j] = 0;
    }

    double Phi[6][6];
    for(int i=0; i<6; i++){
        for(int j=0; j<6; j++){
            Phi[i][j] = 0;
        }
    }

    // Measurement loop
    double t = 0;

    double t_old, Mjd_TT, Mjd_UT1, theta, Azim, Elev, Dist;
    double Y_old[6];
    double Y0[6];
    double U[3][3];
    double K[6][1];
    double x2[6][1]; 
    double P2[6][6];
    double r[3];
    double s[3];
    double Ur[3];
    double UrRs[3];
    double dAds[3];
    double dEds[3];
    double dAdsLT[3];
    double dAdsLTU[3];
    double dAdY[1][6];
    double dEdsLT[3];
    double dEdsLTU[3];
    double dEdY[1][6];
    double dDds[3];
    double dDdsLT[3];
    double dDdsLTU[3];
    double dDdY[1][6];

    for (int i=0; i<nobs; i++){
         // Previous step
        t_old = t;
        Y_old[0] = Y[0][0];
        Y_old[1] = Y[1][0];
        Y_old[2] = Y[2][0];
        Y_old[3] = Y[3][0];
        Y_old[4] = Y[4][0];
        Y_old[5] = Y[5][0];
    
        // Time increment and propagation
        Mjd_UTC = obs[i][0];                       // Modified Julian Date
        t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
    
        IERS (eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;
        
        for (int ii=0; ii<6; ii++){
            yPhi[ii] = Y_old[ii];
            for (int j=0; j<6; j++){
                if (ii == j){
                    yPhi[6*j+ii+6] = 1;
                }else{
                    yPhi[6*j+ii+6] = 0;
                }
            } 
        }
    
        DEIntegVarEqn(0, t-t_old, 1e-13, 1e-6, 42, yPhi);
    
        // Extract state transition matrices
        for (int j = 0; j < 6; j++) {
            for (int i = 0; i < 6; i++) {
                Phi[i][j] = yPhi[6*j + i + 1];
            }
        }
    
        DEIntegAccel(0, t-t_old, 1e-13, 1e-6, 6, Y_old);
        Y[0][0] = Y_old[0];
        Y[1][0] = Y_old[1];
        Y[2][0] = Y_old[2];
        Y[3][0] = Y_old[3];
        Y[4][0] = Y_old[4];
        Y[5][0] = Y_old[5];
    
        // Topocentric coordinates
        theta = gmst(Mjd_UT1);                    // Earth rotation
        R_z(theta, U);
        r[0] = Y[0][0];
        r[1] = Y[1][0];
        r[2] = Y[2][0];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Ur[i] += U[i][j] * r[j];
            }
        }

        UrRs[0] = Ur[0]-Rs[0];
        UrRs[1] = Ur[1]-Rs[1];
        UrRs[2] = Ur[2]-Rs[2];

        // Topocentric position [m]
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                s[i] += LT[i][j] * UrRs[j];
            }
        }

        // Time update
        TimeUpdate(P, Phi);
        
        // Azimuth and partials
        AzElPa(s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dAdsLT[i] += dAds[j] * LT[j][i];
            }
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                dAdsLTU[i] += dAdsLT[j] * U[j][i];
            }
        }

        dAdY[0][0] = dAdsLTU[0];
        dAdY[0][1] = dAdsLTU[1];
        dAdY[0][2] = dAdsLTU[2];
        dAdY[0][3] = 0;
        dAdY[0][4] = 0;
        dAdY[0][5] = 0;

        // Measurement update
        MeasUpdate ( Y, obs[i][1], Azim, sigma_az, dAdY, P, 6, K, x2, P2 ); 
    
        // Elevation and partials
        r[0] = x2[0][0];
        r[1] = x2[1][0];
        r[2] = x2[2][0];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Ur[i] = U[i][j] * r[j];
            }
        }

        UrRs[0] = Ur[0]-Rs[0];
        UrRs[1] = Ur[1]-Rs[1];
        UrRs[2] = Ur[2]-Rs[2];

        // Topocentric position [m]
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                s[i] += LT[i][j] * UrRs[j];
            }
        }

        AzElPa(s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dEdsLT[i] += dEds[j] * LT[j][i];
            }
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                dEdsLTU[i] += dEdsLT[j] * U[j][i];
            }
        }

        dEdY[0][0] = dEdsLTU[0];
        dEdY[0][1] = dEdsLTU[1];
        dEdY[0][2] = dEdsLTU[2];
        dEdY[0][3] = 0;
        dEdY[0][4] = 0;
        dEdY[0][5] = 0;
    
        Y[0][0] = x2[0][0];
        Y[1][0] = x2[1][0];
        Y[2][0] = x2[2][0];
        Y[3][0] = x2[3][0];
        Y[4][0] = x2[4][0];
        Y[5][0] = x2[5][0];
        // Measurement update
        MeasUpdate ( Y, obs[i][3], Elev, sigma_el, dEdY, P, 6, K, x2, P2 );
    
        // Range and partials
        r[0] = x2[0][0];
        r[1] = x2[1][0];
        r[2] = x2[2][0];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Ur[i] = U[i][j] * r[j];
            }
        }

        UrRs[0] = Ur[0]-Rs[0];
        UrRs[1] = Ur[1]-Rs[1];
        UrRs[2] = Ur[2]-Rs[2];

        // Topocentric position [m]
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                s[i] += LT[i][j] * UrRs[j];
            }
        }

        Dist = norm(s, 3); 
        dDds[0] = (s[0]/Dist);         // Range
        dDds[1] = (s[1]/Dist);
        dDds[2] = (s[2]/Dist);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dDdsLT[i] += dDds[j] * LT[j][i];
            }
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                dDdsLTU[i] += dDdsLT[j] * U[j][i];
            }
        }

        dDdY[0][0] = dDdsLTU[0];
        dDdY[0][1] = dDdsLTU[1];
        dDdY[0][2] = dDdsLTU[2];
        dDdY[0][3] = 0;
        dDdY[0][4] = 0;
        dDdY[0][5] = 0;
    
        Y[0][0] = x2[0][0];
        Y[1][0] = x2[1][0];
        Y[2][0] = x2[2][0];
        Y[3][0] = x2[3][0];
        Y[4][0] = x2[4][0];
        Y[5][0] = x2[5][0];
        // Measurement update
        MeasUpdate( Y, obs[i][4], Dist, sigma_range, dDdY, P, 6, K, x2, P2 );
    }  
   

    IERS(eopdata, obs[17][0], UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    double Y6[6];
    Y6[0] = x2[0][0];
    Y6[1] = x2[1][0];
    Y6[2] = x2[2][0];
    Y6[3] = x2[3][0];
    Y6[4] = x2[4][0];
    Y6[5] = x2[5][0];
    DEIntegAccel(0, -(obs[17][0]-obs[0][0])*86400.0, 1e-13, 1e-6, 6, Y6);

    Y0[0] = Y6[0];
    Y0[1] = Y6[1];
    Y0[2] = Y6[2];
    Y0[3] = Y6[3];
    Y0[4] = Y6[4];
    Y0[5] = Y6[5];

    double Y_true[6] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};

    cout<<"Error of Position Estimation"<<endl;

    cout<<"dX%10.1f [m]";
    cout<<Y0[0]-Y_true[0]<<endl;

    cout<<"dY%10.1f [m]";
    cout<<Y0[1]-Y_true[1]<<endl;

    cout<<"dZ%10.1f [m]";
    cout<<Y0[2]-Y_true[2]<<endl;

    cout<<"Error of Velocity Estimation"<<endl;

    cout<<"dVx%8.1f [m/s]";
    cout<<Y0[3]-Y_true[3]<<endl;

    cout<<"dVy%8.1f [m/s]";
    cout<<Y0[4]-Y_true[4]<<endl;

    cout<<"dVz%8.1f [m/s]";
    cout<<Y0[5]-Y_true[5]<<endl;

    return 0;
}
