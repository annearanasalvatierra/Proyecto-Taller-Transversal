#ifndef _MATRIZ_
#define _MATRIZ_

void Transpuesta(double matriz[3][3], double transpuesta[3][3]);
void Producto3X3y3X3(double matriz1[3][3], double matriz2[3][3], double producto[3][3]);
void Producto3X3y3X1(double matriz1[3][3], double matriz2[3][1], double producto[3][1]);
void TripleProducto(double matriz1[3][3], double matriz2[3][3], double matriz3[3][3], double producto[3][3]);

#endif