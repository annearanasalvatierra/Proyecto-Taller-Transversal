/*--------------------------------------------------------------------------

 Purpose:
  Computes azimuth, elevation and partials from local tangent coordinates

 Input:
   s      Topocentric local tangent coordinates (East-North-Zenith frame)
 
 Outputs:
   A      Azimuth [rad]
   E      Elevation [rad]
   dAds   Partials of azimuth w.r.t. s
   dEds   Partials of elevation w.r.t. s

 Last modified:   2015/08/12   M. Mahooti

--------------------------------------------------------------------------*/
#include <cmath>
#include "AzElPa.h"
#include "vector.h"
#include "Sat_Const.h"

void AzElPa(double s[3], double& Az, double& El, double dAds[3], double dEds[3]){

    double pi2 = 2.0 * pi;

    double rho = sqrt(s[0]*s[0] + s[1]*s[1]);

    // Angles
    Az = atan2(s[0],s[1]);

    if (Az < 0.0){
        Az += pi2;
    }

    El = atan ( s[2] / rho );

    // Partials
    dAds[0] = s[1]/(rho*rho);
    dAds[1] = -s[0]/(rho*rho);
    dAds[2] = 0.0 ;

    double dotS = dot(s, s);
    dEds[0] = -s[0]*s[2]/(rho*dotS);
    dEds[1] = -s[1]*s[2]/(rho*dotS);
    dEds[2] = rho/dotS;
}