/*--------------------------------------------------------------------------

 PoleMatrix.m

 Purpose:
   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
   for a given date

 Input:
   Pole coordinte(xp,yp)

 Output:
   PoleMat   Pole matrix

 Last modified:   2015/08/12   M. Mahooti
 
--------------------------------------------------------------------------*/
#include "PoleMatrix.h"
#include "R_x.h"
#include "R_y.h"
#include "Matriz.h"

void PoleMatrix (double xp, double yp, double PoleMat[3][3]){

	double matriz1[3][3] = {0};
	R_y(-xp, matriz1);

	double matriz2[3][3] = {0};
	R_x(-yp, matriz2);

	Producto3X3y3X3(matriz1, matriz2, PoleMat);
}