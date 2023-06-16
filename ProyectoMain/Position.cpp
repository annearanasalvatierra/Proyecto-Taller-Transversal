/*--------------------------------------------------------------------------

 Position.m

 Purpose:
   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
   latitude [rad], altitude [m])

 Last modified:   2015/08/12   M. Mahooti

--------------------------------------------------------------------------*/
#include <cmath>
#include "Sat_Const.h"
#include "Position.h"

void Position(double lon, double lat, double h, double r[3]){
	
	double R_equ = R_Earth;
	double f     = f_Earth;

	double e2     = f*(2.0-f);    // Square of eccentricity
	double CosLat = cos(lat); // (Co)sine of geodetic latitude
	double SinLat = sin(lat);

	// Position vector 
	double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

	r[0] =  (         N+h)*CosLat*cos(lon);
	r[1] =  (         N+h)*CosLat*sin(lon);
	r[2] =  ((1.0-e2)*N+h)*SinLat;
}