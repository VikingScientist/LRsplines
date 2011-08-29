#include <stdio.h>
#include <iostream>
#include <GoTools/geometry/SplineSurface.h>
#include "LRSplineSurface.h"

using namespace Go;
using namespace LR;
using namespace std;

int main() {

	const double TOL = 1e-6;
	
	int p1 = 3;
	int p2 = 3;
	int n1 = 7;
	int n2 = 7;
	int dim = 4;
	bool rat = false;
	double knot_u[n1+p1];
	double knot_v[n2+p2];
	double cp[(dim+rat)*n1*n2];

	// make a uniform integer knot vector
	for(int i=0; i<p1+n1; i++)
		knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
	for(int i=0; i<p2+n2; i++)
		knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;

	// create a list of random control points (all between 0.1 and 1.1)
	int k=0;
	for(int j=0; j<n2; j++) 
		for(int i=0; i<n1; i++) 
			for(int d=0; d<dim+rat; d++)
				cp[k++] = ((i*2+j*3+d*5) % 13) / 13.0 + 0.1;  // rational weights also random and thus we need >0
		
	// make the LR B-spline surface
	LRSplineSurface lr(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);

	// evaluate function values on edges, knots and in between the knots
	for(double u=0; u<=n1-p1+1; u+=0.5) {
		for(double v=0; v<=n2-p2+1; v+=0.5) {

			BasisDerivsSf lr_basis;
			lr.computeBasis(u,v, lr_basis);
			double sum        = 0;
			double sum_diff_u = 0;
			double sum_diff_v = 0;
			for(uint i=0; i<lr_basis.basisValues.size(); i++) {
				sum        += lr_basis.basisValues[i] ;
				sum_diff_u += lr_basis.basisDerivs_u[i];
				sum_diff_v += lr_basis.basisDerivs_v[i];
			}
			cout << "sum (" << u << ", " << v <<") minus one = " << sum-1.0 << endl;
			cout << "sum diff u = " << sum_diff_u << endl;
			cout << "sum diff v = " << sum_diff_v << endl;
			if(fabs(sum-1.0) > TOL || fabs(sum_diff_u) > TOL || fabs(sum_diff_v) > TOL) {
				cout << "ASSERTION FAILED at (" << u << ", " << v << ")\n";
				cout << "=== DIFF U ===\n";
				for(uint i=0; i<lr_basis.basisValues.size(); i++)  {
					if(lr_basis.basisDerivs_u[i] != 0) {
						if(i<10) cout << " ";
						cout << i << ": " << lr_basis.basisDerivs_u[i] << endl;
					}
				}
				cout << "=== DIFF V ===\n";
				for(uint i=0; i<lr_basis.basisValues.size(); i++)  {
					if(lr_basis.basisDerivs_v[i] != 0) {
						if(i<10) cout << " ";
						cout << i << ": " << lr_basis.basisDerivs_v[i] << endl;
					}
				}
			} else
				cout << "ASSERTION PASSED at (" << u << ", " << v << ")\n";
		}
	}
	cout << lr << endl;

}

