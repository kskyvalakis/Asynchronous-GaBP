#ifndef GBP_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define GBP_H

#define N 3		// Declare matrix dimension
#define maxiter 200

double metric(double C[maxiter][N], int r);

double col_sum(int c, double mat[N][N]);

void asynch_GBP(double A[N][N], double b[N], double epsilon, double h[N], double J[N]);

#endif
