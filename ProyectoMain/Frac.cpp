/*--------------------------------------------------------------------------
 
  Fractional part of a number (y=x-[x])

  Last modified:   2015/08/12   M. Mahooti

--------------------------------------------------------------------------*/
#include <cmath>
#include "Frac.h"

double Frac(double x) {
    double res = x - floor(x);
    return res;
}
