#ifndef _GLOBALES_
#define _GLOBALES_

extern double **Snm, **Cnm, **eopdata;

typedef struct{
	double Mjd_TT;
	int n;
	int m;
	double Mjd_UTC;
} Param;

extern Param AuxParam;

#endif