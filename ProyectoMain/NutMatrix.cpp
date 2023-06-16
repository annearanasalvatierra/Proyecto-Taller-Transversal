/*--------------------------------------------------------------------------

 NutMatrix.m

 Purpose:
   Transformation from mean to true equator and equinox

 Input:
   Mjd_TT    Modified Julian Date (Terrestrial Time)

Output:
   NutMat    Nutation matrix

 Last modified:   2015/08/12   M. Mahooti
 
--------------------------------------------------------------------------*/
#include "NutMatrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "Matriz.h"
#include "R_x.h"
#include "R_z.h"

void NutMatrix (double Mjd_TT, double NutMat[3][3]){
	
	// Mean obliquity of the ecliptic
	double eps = MeanObliquity (Mjd_TT);

	// Nutation in longitude and obliquity
	double dpsi, deps;
    NutAngles (Mjd_TT, dpsi, deps);

	// Transformation from mean to true equator and equinox
	double matriz1[3][3]={0};
	R_x(-eps-deps, matriz1);

	double matriz2[3][3]={0};
	R_z(-dpsi, matriz2);

	double matriz3[3][3]={0};
	R_x(+eps, matriz3);

	TripleProducto(matriz1, matriz2, matriz3, NutMat);
}