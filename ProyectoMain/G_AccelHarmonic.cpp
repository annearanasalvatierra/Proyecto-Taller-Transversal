/*--------------------------------------------------------------------------

 G_AccelHarmonic.m

 Purpose:
   Computes the gradient of the Earth's harmonic gravity field 

 Input:
   Mjd_UT      Modified Julian Date (Universal Time)
   r           Satellite position vector in the true-of-date system
   n,m         Gravity model degree and order

 Output:
   G    		Gradient (G=da/dr) in the true-of-date system

 Last modified:   2015/08/12   M. Mahooti

--------------------------------------------------------------------------*/
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"

void G_AccelHarmonic( double r[3], double U[3][3], int n_max, int m_max, double G[3][3]){
    
    double d = 1.0;                // Position increment [m]

    // Iniciamos G con ceros
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j){
            G[i][j] = 0;
        }
        
    }

    double dr[3] = {0,0,0};        //Matriz 3x1 de ceros = Vector de tamaño 3 de ceros

    double da1[3], da2[3], da[3], rdr1[3], rdr2[3];

    // Gradient
    for (int i=0; i<3; i++){
        // Set offset in i-th component of the position vector
        dr[0] = 0;
        dr[1] = 0;
        dr[2] = 0;
        dr[i] = d;
        
        rdr1[0] = r[0] + dr[0]/2;
        rdr1[1] = r[1] + dr[1]/2;
        rdr1[2] = r[2] + dr[2]/2;

        rdr2[0] = r[0] - dr[0]/2;
        rdr2[1] = r[1] - dr[1]/2;
        rdr2[2] = r[2] - dr[2]/2;

        // Acceleration difference
        AccelHarmonic ( rdr1, U, n_max, m_max, da1 );
        AccelHarmonic ( rdr2, U, n_max, m_max, da2 );

        da[0] = da1[0] - da2[0];
        da[1] = da1[1] - da2[1];
        da[2] = da1[2] - da2[2];

        // Derivative with respect to i-th axis
        // Asigna el valor da/d a la columna i de la matriz G 
        G[0][i] = da[0]/d;   
        G[1][i] = da[1]/d;
        G[2][i] = da[2]/d; 
    }
}