#include "MeasUpdate.h"

void MeasUpdate(double x[6][1], double z, double g, double s, double G[1][6], double P[6][6], int n, double K[6][1], double x2[6][1], double P2[6][6]){ 

	
	double Inv_W = s*s;

	double GP[1][6];
	//Calcula el producto de G y P
	for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            GP[0][i] += G[0][j] * P[j][i];
        }
    }

	double transpuestaG[6][1];
	// Calcula la traspuesta de la matriz G
    for (int m = 0; m < 6; m++){
        transpuestaG[m][0] = G[0][m];   
    }

	double GPtranspuestaG;
	//Calculamos el producto de GP por la transpuestaG
	for(int a = 0; a < 6; a++){
			GPtranspuestaG += GP[0][a]*transpuestaG[a][0]; 
	}

	double resultado = 1/(GPtranspuestaG + Inv_W);

	double PTranspuestaG[6][1];
	//Calculamos el producto de P por la transpuestaG
	for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            PTranspuestaG[i][0] += P[i][j] * transpuestaG[j][0];
        }
    }

	for (int j = 0; j < 6; j++) {
            K[j][0] += PTranspuestaG[j][0] * resultado;
    }

	// State update
	double aux[6][1];
	for (int jj = 0; jj < 6; jj++){
		aux[jj][0] = K[jj][0]*(z-g);    
	}

	for (int ii = 0; ii < 6; ii++){
		x2[ii][0] = x[ii][0] + aux[ii][0];    
	}

	// Covariance update

	double aux2[6][6] = {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0},{0, 0, 0, 0, 0, 1}};

	//Multiplicamos K por G
	double KG[6][6];
	for(int t = 0; t < 6; t++){
		for(int l = 0; l < 6; l++){
			KG[t][l] = K[t][0] * G[0][l];
		}
	}

	double resultado2[6][6];
	for(int t = 0; t < 6; t++){
		for(int l = 0; l < 6; l++){
			resultado2[t][l] = aux2[t][l] - KG[l][t];
		}
	}

	//Multiplicamos las matrices resultado2 y P
	for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 6; k++) {
                P2[i][j] += resultado2[i][k] * P[k][j];
            }
        }
    }
}