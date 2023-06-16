#include <cmath>
#include "Matriz.h"

//TRANSPUESTA:

void Transpuesta(double matriz[3][3], double transpuesta[3][3]){
    
    // Calcula la traspuesta de la matriz
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            transpuesta[j][i] = matriz[i][j];
        }    
    }
}

//PRODUCTO DE MATRICES 3X3
     
void Producto3X3y3X3(double matriz1[3][3], double matriz2[3][3], double producto[3][3]){
    
    // Calcula el producto de las dos matrices
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                producto[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }
}
   
//PRODUCTO DE UNA MATRIZ 3X3 y OTRA 3X1

void Producto3X3y3X1(double matriz1[3][3], double matriz2[3][1], double producto[3][1]){
    
    // Realizar el producto de las dos matrices
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 1; j++) {
            for (int k = 0; k < 3; k++) {
                producto[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }
}
   
    //PRODUCTO DE TRES MATRICES 3X3

void TripleProducto(double matriz1[3][3], double matriz2[3][3], double matriz3[3][3], double producto[3][3]){
    
    double producto2[3][3] = {0};

    // Calcula el producto de las dos primeras matrices
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                producto2[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }

    // Calcula el producto de las dos siguientes matrices
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                producto[i][j] += producto2[i][k] * matriz3[k][j];
            }
        }
    }
}