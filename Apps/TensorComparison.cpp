
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <GoTools/geometry/SplineSurface.h>
#include "LRSplineSurface.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {

	// set default parameter values
	const double TOL = 1e-6;
	int p1 = 3;
	int p2 = 4;
	int n1 = 8;
	int n2 = 5;
	int dim = 4;
	bool rat = false;
	string parameters(" parameters: \n" \
	                  "   -p1  <n> polynomial ORDER (degree+1) in first parametric direction\n" \
	                  "   -p2  <n> polynomial order in second parametric direction\n" \
	                  "   -n1  <n> number of basis functions in first parametric direction\n" \
	                  "   -n2  <n> number of basis functions in second parametric direction\n" \
	                  "   -dim <n> dimension of the controlpoints\n" \
	                  "   -help    display (this) help information\n");
	
	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p1") == 0)
			p1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p2") == 0)
			p2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n1") == 0)
			n1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n2") == 0)
			n2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-dim") == 0)
			dim = atoi(argv[++i]);
		else if(strcmp(argv[i], "-help") == 0) {
			cout << "usage: " << argv[0] << endl << parameters;
			exit(0);
		} else {
			cout << "usage: " << argv[0] << endl << parameters;
			exit(1);
		}
	}

	// do some error testing on input
	if(n1 < p1) {
		cerr << "ERROR: n1 must be greater or equal to p1\n";
		exit(2);
	} else if(n2 < p2) {
		cerr << "ERROR: n2 must be greater or equal to p2\n";
		exit(2);
	}


	// make a uniform integer knot vector
	double knot_u[n1+p1];
	double knot_v[n2+p2];
	double cp[(dim+rat)*n1*n2];
	for(int i=0; i<p1+n1; i++)
		knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
	for(int i=0; i<p2+n2; i++)
		knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;

	// create a list of "random" control points (all between 0.1 and 1.1)
	int k=0;
	for(int j=0; j<n2; j++) 
		for(int i=0; i<n1; i++) 
			for(int d=0; d<dim+rat; d++)
				cp[k++] = ((i*2+j*3+d*5) % 13) / 13.0 + 0.1;  // rational weights also random and thus we need >0
		
	// make two identical surfaces
	SplineSurface   ss(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);
	LRSplineSurface lr(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);

	// compare function values on edges, knots and in between the knots
	// as well as all derivatives (up to first derivatives)
	vector<Point> lr_pts(3), ss_pts(3);
	vector<bool> assertion_passed;
	vector<double> par_u_values;
	vector<double> par_v_values;
	for(double u=0; u<=n1-p1+1; u+=0.5) {
		for(double v=0; v<=n2-p2+1; v+=0.5) {
			lr.point(lr_pts, u,v, 1);
			ss.point(ss_pts, u,v, 1, false, false);

			bool correct = true;
			for(int i=0; i<3; i++)  
				for(int d=0; d<dim; d++)
					if( fabs((lr_pts[i][d]-ss_pts[i][d])/ss_pts[i][d]) > TOL ) // relative error
						correct = false;
			cout << "   LR(" << u << ", " << v << ") = " << lr_pts[0] << endl;
			cout << "   SS(" << u << ", " << v << ") = " << ss_pts[0] << endl;
			cout << "dx LR(" << u << ", " << v << ") = " << lr_pts[1] << endl;
			cout << "dx SS(" << u << ", " << v << ") = " << ss_pts[1] << endl;
			cout << "dy LR(" << u << ", " << v << ") = " << lr_pts[2] << endl;
			cout << "dy SS(" << u << ", " << v << ") = " << ss_pts[2] << endl;

			// collect results for summary at the end
			assertion_passed.push_back(correct);
			par_u_values.push_back(u);
			par_v_values.push_back(v);

			if(correct)
				cout << "Parameter (" << u << ", " << v << ") evaluated OK\n";
			else  {
				cout << "ASSERTION FAILED at (" << u << ", " << v << ")\n";
			}
			cout << endl;
		}
	}

	// display results
	bool oneFail = false;
	cout << "     =======    RESULT SUMMARY   ========     \n\n";
	cout << "_____u____________v_____________ASSERT__\n";
	for(uint i=0; i<assertion_passed.size(); i++) {
		printf("%10.4g   %10.4g           ", par_u_values[i], par_v_values[i]);
		if(assertion_passed[i])
			cout << "OK\n";
		else
			cout << "FAIL\n";
		if(!assertion_passed[i])
			oneFail = true;
	}
	cout << "----------------------------------------\n";
	if(oneFail)
		cout << "    test FAILED\n";
	else
		cout << "    all assertions passed\n";
	cout << "========================================\n";

}
