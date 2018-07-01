#include <iostream>
#include <cmath>
#include "gbp.h"

using namespace std;

/* This is an implementation of the (scalar) Gaussian BP algorithm Asynchronous version.
*
************************************************************************************************************
* Input:
*       A - information matrix mxm, (assumed to be symmetric)
*       b - shift vector 1xm
*       max_iter - max iteration
*       epsilon - convergence threshold
*       d - required block size
*
************************************************************************************************************
* Output:
*        The solution for the inference problem
*        vector h of size 1xm s.t. h = max(exp(-1/2x'Ax +x'b)). Equivalently, h = inv(A) * b;
*        J - vector of the precision values Pii (the diagonal of the matrix A^-1)
*        r - round number
*
***************************************************************************************************/
double metric(double C[maxiter][N], int r){
    double norm_num = 0;
    double norm_den = 0;
    double diff[N] = {0};

    for(int i=0; i<N; i++){
        norm_num += pow(C[r][i] - C[r-1][i],2);
        norm_den += pow(C[r][i], 2);
    }
    double div = norm_num/norm_den;

    return div;
}

double col_sum(int c, double mat[N][N]){
    double s = 0;
    for(int i=0; i<N; i++){
        s += mat[i][c];
    }
    return s;
}

void asynch_GBP(double A[N][N], double b[N], double epsilon, double h[N], double J[N]){

    double C[maxiter][N] = {{0}};

    // Initialize messages
    double Mh[N][N] = {{0}};
    double MJ[N][N] = {{0}};

    bool conv = 0;
    double h_j, J_j;

    for(int r=0; r<maxiter; r++){
        for(int i=0; i<N; i++){
            h[i] = b[i] + col_sum(i,Mh);
            J[i] = A[i][i] + col_sum(i,MJ);

            for(int j=0; j<N; j++){
                if(i!=j && A[i][j]!=0){
                    h_j = h[i] - Mh[j][i];
                    J_j = J[i] - MJ[j][i];

                    Mh[i][j] = (-A[j][i]/J_j)*h_j;
                    MJ[i][j] = (-A[j][i]/J_j)*A[i][j];
                }
            }
        }

        for(int i=0; i<N; i++){
            h[i] = b[i] + col_sum(i,Mh);
            J[i] = A[i][i] + col_sum(i,MJ);
        }

        for(int i=0; i<N; i++){
            J[i] = 1/J[i];
            h[i] = h[i]*J[i];
            C[r][i] = h[i];
        }

        if(r>1 && metric(C,r)<epsilon){
            for(int i=0; i<N; i++){
                h[i] = b[i] + col_sum(i,Mh);
                J[i] = A[i][i] + col_sum(i,MJ);
            }
            conv = 1;
            break;
        }
    }

    if(!conv){
        cout << "GaBP did not converge in " << maxiter << " rounds.\n";
        for(int i=0; i<N; i++){
            h[i] = b[i] + col_sum(i,Mh);
            J[i] = A[i][i] + col_sum(i,MJ);
        }
    }

    for(int i=0; i<N; i++){
        J[i] = 1/J[i];
        h[i] = h[i]*J[i];
    }    
}