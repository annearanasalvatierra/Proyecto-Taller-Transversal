#include "Matriz.h"
#include "TimeUpdate.h"

void TimeUpdate(double P[6][6], double Phi[6][6]){
	
	double transpuesta[6][6]={0};
	// Calcula la traspuesta de la matriz
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            transpuesta[j][i] = Phi[i][j];
        }    
    }

	double producto[6][6]={0};
	// Calcula el producto de las dos primeras matrices
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                producto[i][j] += Phi[i][k] * P[k][j];
            }
        }
    }

    // Calcula el producto de las dos siguientes matrices
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                P[i][j] += producto[i][k] * transpuesta[k][j];
            }
        }
    }
}