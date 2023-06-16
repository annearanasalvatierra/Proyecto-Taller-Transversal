/*--------------------------------------------------------------------------

 GHAMatrix.m

 Purpose:
   Transformation from true equator and equinox to Earth equator and 
   Greenwich meridian system 

 Input:
   Mjd_UT1   Modified Julian Date UT1
 
 Output:
   GHAmat    Greenwich Hour Angle matrix

 Last modified:   2015/08/12   M. Mahooti
 
--------------------------------------------------------------------------*/
#include "GHAMatrix.h"
#include "R_z.h"
#include "gast.h"

void GHAMatrix (double Mjd_UT1, double GHAmat[3][3]){

	R_z(gast(Mjd_UT1), GHAmat);
}