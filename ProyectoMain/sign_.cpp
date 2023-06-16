// sign_: returns absolute value of a with sign of b

#include <cmath>
#include "sign_.h"

double sign_(double a, double b){
    
	if (b >= 0.0){
        return fabs(a);
    }else{
        return -fabs(a);
    }
}

