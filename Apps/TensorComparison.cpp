
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <sstream>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/trivariate/SplineVolume.h>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Profiler.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {
#ifdef TIME_LRSPLINE
	Profiler prof(argv[0]);
#endif

	// set default parameter values
	const double TOL = 1e-6;
	int p1 = 3;
	int p2 = 4;
	int p3 = 5;
	int n1 = 8;
	int n2 = 5;
	int n3 = 7;
	int dim = 4;
	bool rat = false;
	bool vol = false;
	stringstream parameters;
	parameters << " parameters: \n" \
	              "   -p1  <n> polynomial ORDER (degree+1) in first parametric direction\n" \
	              "   -p2  <n> polynomial order in second parametric direction\n" \
	              "   -p3  <n> polynomial order in third parametric direction (volume only)\n" \
	              "   -n1  <n> number of basis functions in first parametric direction\n" \
	              "   -n2  <n> number of basis functions in second parametric direction\n" \
	              "   -n2  <n> number of basis functions in third parametric direction (volume only)\n" \
	              "   -dim <n> dimension of the controlpoints\n" \
	              "   -vol     test trivariate volumes instead\n" \
	              "   -help    display (this) help information\n" ;
	parameters << " default values\n";
	parameters << "   -p   = { " << p1 << ", " << p2 << ", " << p3 << " }\n";
	parameters << "   -n   = { " << n1 << ", " << n2 << ", " << n3 << " }\n";
	parameters << "   -vol = { " << ((vol)?"true":"false") << " }\n";

	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p1") == 0)
			p1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p2") == 0)
			p2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p3") == 0) {
			p3  = atoi(argv[++i]);
			vol = true;
		} else if(strcmp(argv[i], "-n1") == 0)
			n1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n2") == 0)
			n2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n3") == 0) {
			n3  = atoi(argv[++i]);
			vol = true;
		} else if(strcmp(argv[i], "-dim") == 0)
			dim = atoi(argv[++i]);
		else if(strcmp(argv[i], "-vol") == 0)
			vol = true;
		else if(strcmp(argv[i], "-help") == 0) {
			cout << "usage: " << argv[0] << endl << parameters.str();
			exit(0);
		} else {
			cout << "usage: " << argv[0] << endl << parameters.str();
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
	double knot_w[n3+p3];
	for(int i=0; i<p1+n1; i++)
		knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
	for(int i=0; i<p2+n2; i++)
		knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;
	for(int i=0; i<p3+n3; i++)
		knot_w[i] = (i<p3) ? 0 : (i>n3) ? n3-p3+1 : i-p3+1;

	// create a list of "random" control points (all between 0.1 and 1.1)
	int nCP = (vol) ? n1*n2*n3 : n1*n2;
	nCP    *= (dim+rat);
	double cp[nCP];
	int k=0;
	for(int i=0; i<nCP; i++) // 839 as a generator over Z_853 gives a period of 425. Should suffice
		cp[k++] = (i*839 % 853) / 853.0 + 0.1;  // rational weights also random and thus we need >0


	bool oneFail = false;
	if(!vol) {
		// make two identical surfaces
		SplineSurface   ss(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);
		LRSplineSurface lr(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);

		lr.generateIDs();

		// compare function values on edges, knots and in between the knots
		// as well as all derivatives (up to first derivatives)
		vector<Point> lr_pts(6), ss_pts(6);
		vector<bool> assertion_passed;
		vector<double> par_u_values;
		vector<double> par_v_values;
		for(double u=0; u<=n1-p1+1; u+=0.5) {
			for(double v=0; v<=n2-p2+1; v+=0.5) {
				lr.point(lr_pts, u,v, 2);
				ss.point(ss_pts, u,v, 2);

				bool correct = true;
				for(int i=0; i<6; i++)
					for(int d=0; d<dim; d++)
						if( fabs(lr_pts[i][d]-ss_pts[i][d]) > TOL ) // absolute error
							correct = false;
				cout << "     LR(" << u << ", " << v << ") = " << lr_pts[0] << endl;
				cout << "     SS(" << u << ", " << v << ") = " << ss_pts[0] << endl;
				cout << "dx   LR(" << u << ", " << v << ") = " << lr_pts[1] << endl;
				cout << "dx   SS(" << u << ", " << v << ") = " << ss_pts[1] << endl;
				cout << "dy   LR(" << u << ", " << v << ") = " << lr_pts[2] << endl;
				cout << "dy   SS(" << u << ", " << v << ") = " << ss_pts[2] << endl;
				cout << "d2x  LR(" << u << ", " << v << ") = " << lr_pts[3] << endl;
				cout << "d2x  SS(" << u << ", " << v << ") = " << ss_pts[3] << endl;
				cout << "dxdy LR(" << u << ", " << v << ") = " << lr_pts[4] << endl;
				cout << "dxdy SS(" << u << ", " << v << ") = " << ss_pts[4] << endl;
				cout << "d2y  LR(" << u << ", " << v << ") = " << lr_pts[5] << endl;
				cout << "d2y  SS(" << u << ", " << v << ") = " << ss_pts[5] << endl;

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
		cout << "     =======    RESULT (SURFACE) SUMMARY   ========     \n\n";
		cout << "_____u____________v_____________ASSERT__\n";
		for(size_t i=0; i<assertion_passed.size(); i++) {
			printf("%10.4g   %10.4g           ", par_u_values[i], par_v_values[i]);
			if(assertion_passed[i])
				cout << "OK\n";
			else
				cout << "FAIL\n";
			if(!assertion_passed[i])
				oneFail = true;
		}

	} else {
		// make two identical volumes
		SplineVolume   sv(n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
		LRSplineVolume lr(n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);

		// compare function values on edges, knots and in between the knots
		// as well as all derivatives (up to first derivatives)
		vector<Point> lr_pts(10), sv_pts(10);
		vector<bool> assertion_passed;
		vector<double> par_u_values;
		vector<double> par_v_values;
		vector<double> par_w_values;
		for(double u=0; u<=n1-p1+1; u+=0.5) {
			for(double v=0; v<=n2-p2+1; v+=0.5) {
				for(double w=0; w<=n3-p3+1; w+=0.5) {
					lr.point(lr_pts, u,v,w,2);
					sv.point(sv_pts, u,v,w,2);

					bool correct = true;
					for(int i=0; i<10; i++)
						for(int d=0; d<dim; d++)
							if( fabs(lr_pts[i][d]-sv_pts[i][d]) > TOL ) // absolute error
								correct = false;
					cout << "     LR(" << u << ", " << v << ") = " << lr_pts[0] << endl;
					cout << "     SS(" << u << ", " << v << ") = " << sv_pts[0] << endl;
					cout << "dx   LR(" << u << ", " << v << ") = " << lr_pts[1] << endl;
					cout << "dx   SS(" << u << ", " << v << ") = " << sv_pts[1] << endl;
					cout << "dy   LR(" << u << ", " << v << ") = " << lr_pts[2] << endl;
					cout << "dy   SS(" << u << ", " << v << ") = " << sv_pts[2] << endl;
					cout << "dy   LR(" << u << ", " << v << ") = " << lr_pts[3] << endl;
					cout << "dy   SS(" << u << ", " << v << ") = " << sv_pts[3] << endl;
					cout << "d2x  LR(" << u << ", " << v << ") = " << lr_pts[4] << endl;
					cout << "d2x  SS(" << u << ", " << v << ") = " << sv_pts[4] << endl;
					cout << "dxdy LR(" << u << ", " << v << ") = " << lr_pts[5] << endl;
					cout << "dxdy SS(" << u << ", " << v << ") = " << sv_pts[5] << endl;
					cout << "dxdz LR(" << u << ", " << v << ") = " << lr_pts[6] << endl;
					cout << "dxdz SS(" << u << ", " << v << ") = " << sv_pts[6] << endl;
					cout << "d2y  LR(" << u << ", " << v << ") = " << lr_pts[7] << endl;
					cout << "d2y  SS(" << u << ", " << v << ") = " << sv_pts[7] << endl;
					cout << "dydz LR(" << u << ", " << v << ") = " << lr_pts[8] << endl;
					cout << "dydz SS(" << u << ", " << v << ") = " << sv_pts[8] << endl;
					cout << "d2z  LR(" << u << ", " << v << ") = " << lr_pts[9] << endl;
					cout << "d2z  SS(" << u << ", " << v << ") = " << sv_pts[9] << endl;

					// collect results for summary at the end
					assertion_passed.push_back(correct);
					par_u_values.push_back(u);
					par_v_values.push_back(v);
					par_w_values.push_back(w);

					if(correct)
						cout << "Parameter (" << u << ", " << v << ") evaluated OK\n";
					else  {
						cout << "ASSERTION FAILED at (" << u << ", " << v << ")\n";
					}
					cout << endl;
				}
			}
		}

		// display results
		cout << "     =======    RESULT (VOLUME) SUMMARY   ========     \n\n";
		cout << "_____u____________v____________w_____________ASSERT__\n";
		for(size_t i=0; i<assertion_passed.size(); i++) {
			printf("%10.4g   %10.4g   %10.4g           ", par_u_values[i], par_v_values[i], par_w_values[i]);
			if(assertion_passed[i])
				cout << "OK\n";
			else
				cout << "FAIL\n";
			if(!assertion_passed[i])
				oneFail = true;
		}

	}
	cout << "----------------------------------------\n";
	if(oneFail)
		cout << "    test FAILED\n";
	else
		cout << "    all assertions passed\n";
	cout << "========================================\n";

}
