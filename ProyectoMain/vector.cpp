#include <cmath>
#include "vector.h"

// Función que calcula la norma de un vector
double norm(double v[], int n){
	
	double suma = 0;
	int i;
	
	if(n <= 0)
		throw "Empty vector";
	
	for(i=0; i<n; i++)
		suma += v[i]*v[i];
		
	return(sqrt(suma));
}

// Función que calcula el producto escalar de dos vectores
double dot(double v1[], double v2[], int nv1, int nv2)
{
    double sum = 0;
    int i;

    if((nv1 <= 0) || (nv2 <= 0) || (nv1 != nv2))
        throw "Empty vector";

    for(i = 0; i < nv1; i++)
        sum += v1[i]*v2[i];

    return sum;
}

// Función que calcula el producto vectorial de dos vectores
void cross(double v[], int &n, double v1[], double v2[], int nv1, int nv2)
{
    double sum = 0;

    if((nv1 <= 0) || (nv2 <= 0) || (nv1 != nv2))
        throw "Empty vector";

    v[0] = -v1[2]*v2[1] + v1[1]*v2[2];
    v[1] =  v1[2]*v2[0] - v1[0]*v2[2];
    v[2] = -v1[1]*v2[0] + v1[0]*v2[1];
    n = nv1;
}