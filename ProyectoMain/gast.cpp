/*--------------------------------------------------------------------------

 GAST.m

 Purpose:
   Greenwich Apparent Sidereal Time

 Input:
   Mjd_UT1   Modified Julian Date UT1

 Output:
   gstime    GAST in [rad]
 
 Last modified:   2015/08/12   M. Mahooti
 
--------------------------------------------------------------------------*/
#include <cmath>
#include "gast.h"
#include "gmst.h"
#include "EqnEquinox.h"
#include "Sat_Const.h"

double gast (double Mjd_UT1){

	return fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*pi );
}

