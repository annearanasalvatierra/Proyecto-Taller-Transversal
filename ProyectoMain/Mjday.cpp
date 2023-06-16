/*--------------------------------------------------------------------------
  inputs:
    year        - year                           
    mon         - month                          
    day         - day                            
    hr          - universal time hour            
    min         - universal time min             
    sec         - universal time sec             

  output:       
    Mjd         - Modified julian date                    
--------------------------------------------------------------------------*/

#include <cmath>
#include "Mjday.h"

double Mjday(int yr, int mon, int day, int hr, int min, double sec){

    double jd = 367.0 * yr - floor( (7.0 * (yr + floor( (mon + 9.0) / 12.0) ) ) * 0.25 ) + floor( 275.0 * mon / 9.0 ) + day + 1721013.5 + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;

    return jd-2400000.5;
}