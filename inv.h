#ifndef INV_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define INV_H

#define N 3		// Declare matrix dimension

bool inverse(double A[N][N], double inverse[N][N]);

void getCofactor(double A[N][N], double temp[N][N], int p, int q, int n);

double determinant(double A[N][N]);

void adjoint(double A[N][N],double adj[N][N]);

#endif