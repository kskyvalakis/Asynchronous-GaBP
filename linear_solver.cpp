#include <iostream>
#include "inv.h"
#include "display.h"
#include "gbp.h"

using namespace std;

void multiply(double b[N], double inv[N][N], double m[N])
{
   int i, j, k;
   for (i = 0; i < N; i++){
   		for(j = 0; j < N; j++){
   			m[i] += b[j]*inv[i][j];
   		}
   }
}

int main()
{
	double b[N] = {1, 1, 1};
	double A[N][N] = {{1, 0.3, 0.1},{0.3, 1, 0.1},{0.1, 0.1, 1}};

	// Solve the syste using the convential inverse matrix method
	double inv[N][N], x[N];			// Store inverse matrix and solution vector x = A^-1*b
	inverse(A,inv);
	multiply(b,inv,x);
	cout << "\nConventional solution : ";
	for (int j=0; j<N; j++)
    	cout << x[j] << " ";
    cout << endl;

    cout << "\nConventional solution matrix diagonal : ";
	for (int j=0; j<N; j++)
    	cout << inv[j][j] << " ";
    cout << endl;


    // Solve the system using GaBP algorithm
    double epsilon = 0.000001;
    double h[N] = {0};
    double J[N] = {0};
    asynch_GBP(A, b, epsilon, h, J);

    cout << "\nGaBP solution : ";
	for (int j=0; j<N; j++)
    	cout << h[j] << " ";
    cout << endl;

    cout << "\nGaBP solution matrix diagonal : ";
	for (int j=0; j<N; j++)
    	cout << J[j] << " ";
    cout << endl;

}