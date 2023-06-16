#include "Legendre.h"
#include <cmath>
#include <iostream>

using namespace std;

void Legendre(int n, int m, double fi, double*** pnm, double*** dpnm){

    (*pnm) = new double*[n+1];
    (*dpnm) = new double*[n+1];

    for (int i = 0; i < n+1; i++) { 
        (*pnm)[i]  = new double[m+1];
        (*dpnm)[i] = new double[m+1];
    }

    // Iniciamos a cero
    for (int i = 0; i < n+1; i++) { 
        for(int j=0; j<m+1; j++){
            (*pnm)[i][j]  = 0;
            (*dpnm)[i][j] = 0;
        } 
    }

    (*pnm)[0][0] =1;
    (*dpnm)[0][0]=0;
    (*pnm)[1][1] =sqrt(3)*cos(fi);
    (*dpnm)[1][1]=-sqrt(3)*sin(fi);

    // diagonal coefficients
    for (int i=2; i<n+1; i++){
        (*pnm)[i][i] = sqrt((2.0*i+1.0)/(2.0*i))*cos(fi)*(*pnm)[i-1][i-1];
    }   
    
    for (int j=2; j<n+1; j++){
         (*dpnm)[j][j] = sqrt((2.0*j+1.0)/(2.0*j))*((cos(fi)*(*dpnm)[j-1][j-1])-(sin(fi)*(*pnm)[j-1][j-1]));
    }
       
    // horizontal first step coefficients
    for (int k=1; k<n+1; k++){
        (*pnm)[k][k-1]= sqrt(2*k+1)*sin(fi)*(*pnm)[k-1][k-1];
    }
       
    for (int s=1; s<n+1; s++){
        (*dpnm)[s][s-1]= sqrt(2*s+1)*((cos(fi)*(*pnm)[s-1][s-1])+(sin(fi)*(*dpnm)[s-1][s-1]));
    }
        
    // horizontal second step coefficients
    int i = 0;
    int j = 0;
    int k = 2;
    while(1){
        for (i=k; i<n+1; i++){
            (*pnm)[i][j]=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*(*pnm)[i-1][j])-(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3.0))*(*pnm)[i-2][j]));
        }        
        j = j+1;
        k = k+1;
        if (j>m){
            break;
        }
    }
        
    j = 0;
    k = 2;
    while(1){
        for (i=k; i<n+1; i++){
            (*dpnm)[i][j]=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*(*dpnm)[i-1][j])+(sqrt(2*i-1)*cos(fi)*(*pnm)[i-1][j])-(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3.0))*(*dpnm)[i-2][j]));
        }      
        j = j+1;
        k = k+1;
        if (j>m){
            break;
        }
    }
}