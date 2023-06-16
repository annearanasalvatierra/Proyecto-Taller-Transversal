#include <cmath>
#include <stdexcept>
#include "Geodetic.h"
#include "Sat_Const.h"
#include "vector.h"

void Geodetic( double r[3], double& lon, double& lat, double& h){
	
	double R_equ = R_Earth;
	double f     = f_Earth;

	double epsRequ = 2.220446049250313e-16*R_equ;        // Convergence criterion
	double e2      = f*(2.0-f);        // Square of eccentricity
	
	double X = r[0];                   // Cartesian coordinates
	double Y = r[1];
	double Z = r[2];
	double rho2 = X*X + Y*Y;           // Square of distance from z-axis
	
	// Check validity of input data
	if (norm(r, 3)==0.0){
		throw std::runtime_error("Invalid input in Geodetic constructor");
		lon = 0.0;
		lat = 0.0;
		h   = -R_Earth;
	}
	
	// Iteration 
	double dZ = e2*Z;
	
	double ZdZ, Nh, SinPhi, N, dZ_new;

	while(1){
		ZdZ    =  Z + dZ;
		Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
		SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
		N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
		dZ_new =  N*e2*SinPhi;
		if ( abs(dZ-dZ_new) < epsRequ ){
			break;
		}
		dZ = dZ_new;
	}
	
	// Longitude, latitude, altitude
	lon = atan2 ( Y, X );
	lat = atan2 ( ZdZ, sqrt(rho2) );
	h   = Nh - N;
}