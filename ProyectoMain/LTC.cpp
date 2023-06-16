/*--------------------------------------------------------------------------
 
 LTC.m

 Purpose:
   Transformation from Greenwich meridian system to 
   local tangent coordinates

 Inputs:
   lon      -Geodetic East longitude [rad]
   lat      -Geodetic latitude [rad]
   
 Output:
   M        -Rotation matrix from the Earth equator and Greenwich meridian
             to the local tangent (East-North-Zenith) coordinate system

 Last modified:   2015/08/12   M. Mahooti

--------------------------------------------------------------------------*/
#include "LTC.h"
#include "R_y.h"
#include "R_z.h"
#include "Matriz.h"

void LTC(double lon, double lat, double M[3][3]){
    
    double matriz1[3][3]={0};
    R_y(-1.0*lat, matriz1);

    double matriz2[3][3]={0};
    R_z(lon, matriz2);

    Producto3X3y3X3(matriz1, matriz2, M);
    
    double Aux;
    for(int k=0; k<3; k++){
        Aux = M[0][k];
        M[0][k] = M[1][k];
        M[1][k] = M[2][k];
        M[2][k] = Aux;
    }
}